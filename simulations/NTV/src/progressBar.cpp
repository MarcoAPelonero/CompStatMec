#include "ProgressBar.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

ProgressBar::ProgressBar(int total, int width, const std::string& prefix, const std::string& suffix)
    : totalSteps(total), currentStep(0), barWidth(width), prefixText(prefix), suffixText(suffix) {
    startTime = std::chrono::steady_clock::now();
    display();
}

ProgressBar::~ProgressBar() {
    finish();
}

void ProgressBar::step() {
    if (currentStep < totalSteps) {
        currentStep++;
        display();
    }
}

void ProgressBar::update(int current) {
    if (current >= currentStep && current <= totalSteps) {
        currentStep = current;
        display();
    }
}

void ProgressBar::finish() {
    if (currentStep < totalSteps) {
        currentStep = totalSteps;
        display();
    }
    std::cout << std::endl;
}

void ProgressBar::display() {
    float progress = static_cast<float>(currentStep) / totalSteps;
    int pos = static_cast<int>(barWidth * progress);

    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - startTime).count();

    std::cout << "\r" << prefixText << " [";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "% "
              << suffixText << " Elapsed: " << elapsed << "s";
    std::cout.flush();
}
