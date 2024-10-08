#include "progress_bar.hpp"   // Include the header file
#include <chrono>
#include <thread>
#include <iomanip>             // For std::setprecision

// Constructor Implementation
ProgressBar::ProgressBar(int total_steps, int bar_width, std::string bar_char)
    : total_steps(total_steps), bar_width(bar_width), bar_char(bar_char) {
    start_time = std::chrono::steady_clock::now();
}

// Update Function Implementation
void ProgressBar::update(int step) {
    // Calculate progress
    double progress = static_cast<double>(step) / total_steps;
    int pos = static_cast<int>(bar_width * progress);

    // Calculate time elapsed
    auto current_time = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

    // Display progress bar
    std::cout << "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos)
            std::cout << bar_char;
        else
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << (progress * 100.0) << "% ";
    std::cout << "Elapsed: " << time_elapsed << "s\r";
    std::cout.flush();
}

// Finish Function Implementation
void ProgressBar::finish() {
    update(total_steps);
    std::cout << std::endl;
}
