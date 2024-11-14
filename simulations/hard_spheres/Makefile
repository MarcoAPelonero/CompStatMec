# Compiler
CXX := g++
# Compiler flags
CXXFLAGS := -Wall -Wextra -std=c++17
# Directories for include files and object files
INCLUDES := .
OBJ_DIR := obj

# Source files
SOURCES := main.cpp scalar.cpp vector.cpp
# Object files
OBJECTS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
# Target executable
TARGET := main_executable

# Rule to build the target
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile source files into object files
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c $< -o $@

# Create object directory if it doesn't exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Clean up build files
clean:
	rm -rf $(OBJ_DIR) $(TARGET)
