#ifndef PROGRESS_BAR_HPP
#define PROGRESS_BAR_HPP

#include <iostream>
#include <chrono>
#include <string>

// ProgressBar class definition
class ProgressBar {
private:
    int total_steps;           // Total number of steps in the process
    int bar_width;             // Width of the progress bar
    std::string bar_char;      // Character to represent the progress
    std::chrono::steady_clock::time_point start_time; // Start time for tracking elapsed time

public:
    // Constructor
    ProgressBar(int total_steps, int bar_width = 50, std::string bar_char = "#");

    // Update the progress bar
    void update(int step);

    // Complete the progress bar
    void finish();
};

#endif // PROGRESS_BAR_HPP