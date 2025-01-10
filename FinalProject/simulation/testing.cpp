#include <iostream>
#include <fstream>
#include <vector>
#include "include/json.hpp"

using json = nlohmann::json;

int main() {
    // Create proper JSON structure
    std::vector<float> rdf(1000, 0.5f);
    json test_json = {
        {"rdf", rdf},
        {"density", 0.8f}
    };

    // Write to file
    std::ofstream json_file("test_input.json");
    json_file << test_json.dump(2);
    json_file.close();

    // Execute Python script
    std::string command = "python scripts/predict_energy.py layers/best_model.pth < test_json.json";
    FILE* pipe = _popen(command.c_str(), "r");
    
    if (!pipe) {
        std::cerr << "Error executing command" << std::endl;
        return 1;
    }

    // Read output
    char buffer[128];
    std::string result;
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

    int status = _pclose(pipe);
    
    std::cout << "Exit status: " << status << std::endl;
    std::cout << "Output: " << result << std::endl;

    return 0;
}