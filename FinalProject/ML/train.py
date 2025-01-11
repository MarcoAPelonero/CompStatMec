import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm

from lib.data import Dataset as RDFDatasetLoader
from lib.models import RDFDensityDataset, RDFDensityMLP


def rmse(y_pred, y_true):
    """Root Mean Square Error (RMSE)."""
    return torch.sqrt(torch.mean((y_pred - y_true)**2))


def pick_best_model(trained_models, df_results):
    """
    Selects the single best model with the lowest validation loss.

    Args:
        trained_models (list): List of trained RDFDensityMLP models.
        df_results (pd.DataFrame): DataFrame containing training logs.

    Returns:
        nn.Module: The best-performing model.
    """
    best_val_loss_per_model = df_results.groupby('model_id')['val_loss'].min()

    best_model_id = best_val_loss_per_model.idxmin()
    best_model_index = best_model_id - 1  
    best_model = trained_models[best_model_index]

    return best_model


def train_one_model(model_idx,
                    model,
                    train_loader,
                    val_loader,
                    epochs,
                    optimizer,
                    scheduler,  
                    criterion,
                    device,
                    main_bar=None):
    """
    Train a single model and return epoch-wise logs.
    Use a sub-tqdm bar for this model's epochs.
    After training finishes, that bar is closed and we optionally update the main bar.
    """
    epoch_logs = []  

    with tqdm(total=epochs, desc=f"Training Model {model_idx+1}", leave=False) as pbar:
        best_val_loss = float('inf') 

        for epoch in range(epochs):
            model.train()
            running_loss = 0.0
            running_rmse = 0.0

            for X_rdf, X_dens, y_true in train_loader:
                X_rdf = X_rdf.to(device, dtype=torch.float)       
                X_dens = X_dens.to(device, dtype=torch.float)     
                y_true = y_true.to(device, dtype=torch.float)     

                optimizer.zero_grad()
                y_pred = model(X_rdf, X_dens)
                loss = criterion(y_pred, y_true)
                loss.backward()
                optimizer.step()

                running_loss += loss.item() * X_rdf.size(0)
                running_rmse += rmse(y_pred, y_true).item() * X_rdf.size(0)

            epoch_train_loss = running_loss / len(train_loader.dataset)
            epoch_train_rmse = running_rmse / len(train_loader.dataset)

            model.eval()
            running_val_loss = 0.0
            running_val_rmse = 0.0
            with torch.no_grad():
                for X_rdf_val, X_dens_val, y_val in val_loader:
                    X_rdf_val = X_rdf_val.to(device, dtype=torch.float)
                    X_dens_val = X_dens_val.to(device, dtype=torch.float)
                    y_val = y_val.to(device, dtype=torch.float)

                    y_pred_val = model(X_rdf_val, X_dens_val)
                    val_loss = criterion(y_pred_val, y_val).item()
                    val_rmse_val = rmse(y_pred_val, y_val).item()

                    running_val_loss += val_loss * X_rdf_val.size(0)
                    running_val_rmse += val_rmse_val * X_rdf_val.size(0)

            epoch_val_loss = running_val_loss / len(val_loader.dataset)
            epoch_val_rmse = running_val_rmse / len(val_loader.dataset)

            if epoch_val_loss < best_val_loss:
                best_val_loss = epoch_val_loss
                best_model_state = model.state_dict()

            if scheduler is not None and (epoch + 1) % 10 == 0:
                scheduler.step()

            epoch_logs.append((
                epoch+1,
                epoch_train_loss, epoch_train_rmse,
                epoch_val_loss, epoch_val_rmse
            ))

            pbar.set_postfix({
                "Train_Loss": f"{epoch_train_loss:.4f}",
                "Val_Loss": f"{epoch_val_loss:.4f}"
            })
            pbar.update(1)

    model.load_state_dict(best_model_state)

    if main_bar is not None:
        main_bar.update(1)

    return epoch_logs


def main():
    dataset_folder = "dataset"
    data_obj = RDFDatasetLoader(dataset_folder=dataset_folder)
    data_obj.read_dataset()
    data_obj.split(train=0.8, vali=0.1, test=0.1, random_state=42)

    print(f"Training samples: {len(data_obj.train)}")
    print(f"Validation samples: {len(data_obj.vali)}")
    print(f"Testing samples: {len(data_obj.test)}\n")

    max_len = 1000 
    train_data = RDFDensityDataset(samples=data_obj.train, max_len=max_len)
    val_data   = RDFDensityDataset(samples=data_obj.vali, max_len=max_len)
    test_data  = RDFDensityDataset(samples=data_obj.test, max_len=max_len)

    batch_size = 16
    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    val_loader   = DataLoader(val_data,   batch_size=batch_size, shuffle=False)
    test_loader  = DataLoader(test_data,  batch_size=batch_size, shuffle=False)

    rdf_size = max_len  
    hidden_dim = 128
    lr = 1e-3
    epochs = 80

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}\n")

    results_list = []
    trained_models = []

    overall_best_val_loss = float('inf')
    overall_best_model = None

    try:
        n_models = int(input("Number of models to train: "))
        if n_models < 1:
            raise ValueError
    except ValueError:
        print("Please enter a valid positive integer for the number of models.")
        return
    print(f"Will train {n_models} model(s).")

    with tqdm(total=n_models, desc='Overall Progress', leave=True) as main_bar:
        for m_idx in range(n_models):
            model = RDFDensityMLP(rdf_size=rdf_size, hidden_dim=hidden_dim).to(device)
            criterion = nn.MSELoss()
            optimizer = optim.Adam(model.parameters(), lr=lr)

            scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=10, gamma=0.8, verbose=False)

            epoch_logs = train_one_model(
                model_idx=m_idx,
                model=model,
                train_loader=train_loader,
                val_loader=val_loader,
                epochs=epochs,
                optimizer=optimizer,
                scheduler=scheduler, 
                criterion=criterion,
                device=device,
                main_bar=main_bar
            )

            for (epoch_id, tr_loss, tr_rmse, va_loss, va_rmse) in epoch_logs:
                results_list.append({
                    'model_id': m_idx + 1,
                    'epoch': epoch_id,
                    'train_loss': tr_loss,
                    'train_rmse': tr_rmse,
                    'val_loss': va_loss,
                    'val_rmse': va_rmse
                })

            min_val_loss = min([log[3] for log in epoch_logs]) 
            if min_val_loss < overall_best_val_loss:
                overall_best_val_loss = min_val_loss
                overall_best_model = model

            trained_models.append(model)

    df_results = pd.DataFrame(results_list)

    sns.set(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sns.lineplot(data=df_results, x='epoch', y='train_loss', hue='model_id', ax=axes[0])
    axes[0].set_title("Train Loss vs. Epoch")
    axes[0].set_xlabel("Epoch")
    axes[0].set_ylabel("Train MSE Loss")

    sns.lineplot(data=df_results, x='epoch', y='val_loss', hue='model_id', ax=axes[1])
    axes[1].set_title("Validation Loss vs. Epoch")
    axes[1].set_xlabel("Epoch")
    axes[1].set_ylabel("Validation MSE Loss")

    plt.tight_layout()
    plt.savefig("training_results_seaborn.png", dpi=150)
    plt.show()

    if overall_best_model is not None:
        overall_best_model.eval()
        test_loss = 0.0
        test_rmse_val = 0.0
        with torch.no_grad():
            for X_rdf_test, X_dens_test, y_test in test_loader:
                X_rdf_test = X_rdf_test.to(device, dtype=torch.float)
                X_dens_test = X_dens_test.to(device, dtype=torch.float)
                y_test = y_test.to(device, dtype=torch.float)

                y_pred_test = overall_best_model(X_rdf_test, X_dens_test)
                loss = criterion(y_pred_test, y_test).item()
                test_loss += loss * X_rdf_test.size(0)
                test_rmse_val += rmse(y_pred_test, y_test).item() * X_rdf_test.size(0)

        test_loss /= len(test_loader.dataset)
        test_rmse_val /= len(test_loader.dataset)
        print(f"\n[Best Model Test] MSE: {test_loss:.4f}, RMSE: {test_rmse_val:.4f}")
    else:
        print("Skipping test evaluation as no valid model was trained.")

    if overall_best_model is not None:
        overall_best_model.eval()
        train_loss = 0.0
        train_rmse_val = 0.0
        with torch.no_grad():
            for X_rdf_train, X_dens_train, y_train in train_loader:
                X_rdf_train = X_rdf_train.to(device, dtype=torch.float)
                X_dens_train = X_dens_train.to(device, dtype=torch.float)
                y_train = y_train.to(device, dtype=torch.float)

                y_pred_train = overall_best_model(X_rdf_train, X_dens_train)
                loss = criterion(y_pred_train, y_train).item()
                train_loss += loss * X_rdf_train.size(0)
                train_rmse_val += rmse(y_pred_train, y_train).item() * X_rdf_train.size(0)

        train_loss /= len(train_loader.dataset)
        train_rmse_val /= len(train_loader.dataset)

        val_loss = 0.0
        val_rmse_val = 0.0
        with torch.no_grad():
            for X_rdf_val, X_dens_val, y_val in val_loader:
                X_rdf_val = X_rdf_val.to(device, dtype=torch.float)
                X_dens_val = X_dens_val.to(device, dtype=torch.float)
                y_val = y_val.to(device, dtype=torch.float)

                y_pred_val = overall_best_model(X_rdf_val, X_dens_val)
                loss = criterion(y_pred_val, y_val).item()
                val_loss += loss * X_rdf_val.size(0)
                val_rmse_val += rmse(y_pred_val, y_val).item() * X_rdf_val.size(0)

        val_loss /= len(val_loader.dataset)
        val_rmse_val /= len(val_loader.dataset)

        # Print LaTeX Table
        table = f"""
                \\begin{{table}}[h!]
                    \\centering
                    \\caption{{Comparison of Metrics for the Best Model}}
                    \\label{{tab:metric_comparison}}
                    \\begin{{tabular}}{{lcc}}
                        \\toprule
                        \\textbf{{Dataset}} & \\textbf{{MSE}} & \\textbf{{RMSE}} \\\\
                        \\midrule
                        Train & {train_loss:.4f} & {train_rmse_val:.4f} \\\\
                        Validation & {val_loss:.4f} & {val_rmse_val:.4f} \\\\
                        Test & {test_loss:.4f} & {test_rmse_val:.4f} \\\\
                        \\bottomrule
                    \\end{{tabular}}
                \\end{{table}}
                """
        print(table)
    else:
        print("No metrics to display as no valid model was trained.")

if __name__ == "__main__":
    main()
