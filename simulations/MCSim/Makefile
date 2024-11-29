# Project directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
TEST_DIR = test

# Compiler settings
CC = g++  # Use g++ for compiling both C++ files and C++ test files
CXX = g++  # For compiling C++ files
CFLAGS = -Wall -I$(INCLUDE_DIR)  # -I to include the header directory for C++ files
LDFLAGS =  # Any necessary linker flags can be added here

# Find all source files in the src directory (only .cpp files)
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)  # Object files should be created from .cpp files
EXEC = $(BUILD_DIR)/program

# Find all test files in the test directory (only .cpp files)
TESTS = $(wildcard $(TEST_DIR)/*.cpp)

# Default target
all: $(EXEC)

# Create the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Rule to build the program
$(EXEC): $(OBJS)
	$(CC) $(OBJS) -o $(EXEC) $(LDFLAGS)

# Rule to compile each source file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Run the compiled program
run: $(EXEC)
	$(EXEC)

# Test a specific test file (using TEST variable passed to make)
.PHONY: test
test:
	@if [ -z "$(TEST)" ]; then \
		echo "Error: No test specified. Usage: make test TEST=test_name"; \
		exit 1; \
	fi
	@echo "Running test: $(TEST)"
	# Ensure the object files are up-to-date
	$(MAKE) $(BUILD_DIR)/particles.o
	$(MAKE) $(BUILD_DIR)/vec.o
	$(CXX) $(CFLAGS) $(TEST_DIR)/$(TEST).cpp -o $(BUILD_DIR)/$(TEST) $(BUILD_DIR)/particles.o $(BUILD_DIR)/vec.o && ./$(BUILD_DIR)/$(TEST)

# Run all tests (non-interactive)
.PHONY: test_all
test_all:
	@for test in $(basename $(notdir $(TESTS))); do \
		echo "Running test: $$test"; \
		$(CXX) $(CFLAGS) $(TEST_DIR)/$$test.cpp -o $(BUILD_DIR)/$$test $(BUILD_DIR)/particles.o $(BUILD_DIR)/vec.o && ./$(BUILD_DIR)/$$test; \
	done
