from pathlib import Path
from tqdm import tqdm
import os 
import matplotlib.pyplot as plt
import numpy as np

datasetFolder = "dataset"

class RDF(object):
    def __init__(self, filename):
        self.filename = filename 

    def readFile(self):
        parts = self.filename.split('_')
        rho_index = parts.index('rho') + 1
        T_index = parts.index('T') + 1
        
        r = []
        rdf = []

        self.rho = float(parts[rho_index])
        self.T = float(parts[T_index][:3])

        filepath = os.path.join(datasetFolder, self.filename)
        with open(filepath, 'r') as f:
            for line in f.readlines():
                if line.startswith('# Average Pressure'):
                    self.avg_pressure = float(line.split()[-1])
                    continue
                if line.startswith('# Average Potential'):
                    self.avg_energy = float(line.split()[-1])
                    continue
                if line.startswith('# r	g(r)'):
                    continue
                r.append(float(line.split()[0]))
                rdf.append(float(line.split()[1]))
        self.r = np.array(r) 
        self.rdf = np.array(rdf)

    def show(self, single = True):
        plt.plot(self.r, self.rdf, '-') 
        if single:
            plt.grid(True)
            plt.show()
                
class Dataset(object):
    def __init__(self, datasetfolder, filetype):
        self.filetype = filetype
        self.datasetFolder = datasetfolder
    
    def readDataset(self):
        data = []
        files = list(Path(self.datasetFolder).glob(self.filetype))
        for rdf_file in tqdm(files, desc='Processing RDF Data'):
            filename = rdf_file.stem
            element = RDF(filename)
            element.readFile()
            data.append(element)

        self.data = np.array(data)

    def split(self, train = 0.8, vali = 0.0, test = 0.2):
        if (train+vali+test != 1.0):
            raise ValueError ("Perchantages choosen for the split are not valid")

        dataLenght = len(self.data)
        trainLenght = int(train * dataLenght)
        valiLenght = int(vali * dataLenght)
        testLenght = int(test * dataLenght)

        while (trainLenght + valiLenght + testLenght != dataLenght):
            trainLenght += 1

        for i, dat in enumerate(self.data):
            j  = np.random.randint(0, dataLenght-1)
            self.data[i] = self.data[j]
            self.data[j] = dat
        
        self.train = self.data[:trainLenght]
        self.test = self.data[:testLenght]
        

        
avg_rdf_files = list(Path(datasetFolder).glob('avg_rdf_rho_*_T_*.dat'))
rdf_files = list(Path(datasetFolder).glob('rdf_rho_*_T_*.dat'))

for rdf_file in tqdm(rdf_files, desc='Processing RDF Data'):
    filename = rdf_file.stem
    parts = filename.split('_')
    rho_index = parts.index('rho') + 1
    T_index = parts.index('T') + 1

diocane = RDF("avg_rdf_rho_0.3_T_1.4.dat")
diocane.readFile()