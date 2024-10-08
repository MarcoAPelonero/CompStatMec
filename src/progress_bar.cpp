#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>  // For std::setprecision

class ProgressBar {
private:
    int total_steps;
    int bar_width;
    std::string bar_char;
    std::chrono::steady_clock::time_point start_time;

public:
    ProgressBar(int total_steps, int bar_width = 50, std::string bar_char = "#")
        : total_steps(total_steps), bar_width(bar_width), bar_char(bar_char) {
        start_time = std::chrono::steady_clock::now();
    }

    void update(int step) {
        // Calculate progress
        double progress = (double)step / total_steps;
        int pos = static_cast<int>(bar_width * progress);

        // Calculate time elapsed
        auto current_time = std::chrono::steady_clock::now();
        auto time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

        // Display progress bar
        std::cout << "[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << bar_char;
            else std::cout << " ";
        }
        std::cout << "] " << std::setprecision(2) << std::fixed << (progress * 100.0) << "% ";
        std::cout << "Elapsed: " << time_elapsed << "s\r";
        std::cout.flush();
    }

    void finish() {
        update(total_steps);
        std::cout << std::endl;
    }
};
