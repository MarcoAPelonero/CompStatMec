#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

int main() {
    // Command to invoke the Python script
    const std::string pythonScript = "python dummy_script.py";

    // Open a pipe to the Python script
    FILE* pipe = popen(pythonScript.c_str(), "w");
    if (!pipe) {
        std::cerr << "Failed to open pipe to Python script.\n";
        return 1;
    }

    // Send data to the Python script
    std::string dataToSend = "Hello from C++!";
    fprintf(pipe, "%s\n", dataToSend.c_str());
    fflush(pipe); // Ensure the data is sent immediately

    // Close the pipe
    pclose(pipe);

    std::cout << "Data sent to Python script: " << dataToSend << "\n";

    return 0;
}
