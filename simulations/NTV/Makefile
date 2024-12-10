SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
TEST_DIR = test

# Compiler and flags
CC = g++
CXX = g++
CFLAGS = -g -Wall -I$(INCLUDE_DIR)
EXTRA_CFLAGS = -Wextra -pedantic -Wconversion -Wshadow
LDFLAGS =

# Source files and object files
MAIN_SRC = $(SRC_DIR)/main.cpp
SRCS = $(filter-out $(MAIN_SRC), $(wildcard $(SRC_DIR)/*.cpp))
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
MAIN_OBJ = $(BUILD_DIR)/main.o
EXEC = $(BUILD_DIR)/main  # The target executable

# Test files
TESTS = $(wildcard $(TEST_DIR)/*.cpp)

# Default target
all: $(EXEC)

# Create build directory if it doesn't exist
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Compile the program
$(EXEC): $(OBJS) $(MAIN_OBJ)
	$(CC) $(OBJS) $(MAIN_OBJ) -o $(EXEC) $(LDFLAGS)

# Compile individual object files with default warnings
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile main object file with default warnings
$(MAIN_OBJ): $(MAIN_SRC) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build directory
clean:
	rm -rf $(BUILD_DIR)

# Run the simulation
run: $(EXEC)
	$(EXEC)

# Run a specific test with extra warnings
.PHONY: test
test:
	@if [ -z "$(TEST)" ]; then \
		echo "Error: No test specified. Usage: make test TEST=test_name"; \
		exit 1; \
	fi
	@echo "Running test: $(TEST)"
	# Ensure all object files are up-to-date
	$(MAKE) $(OBJS)
	$(CXX) $(CFLAGS) $(EXTRA_CFLAGS) -c $(TEST_DIR)/$(TEST).cpp -o $(BUILD_DIR)/$(TEST).o
	$(CXX) $(BUILD_DIR)/$(TEST).o $(OBJS) -o $(BUILD_DIR)/$(TEST) $(LDFLAGS) && ./$(BUILD_DIR)/$(TEST)

# Run all tests with extra warnings
.PHONY: test_all
test_all:
	@for test in $(basename $(notdir $(TESTS))); do \
		echo "Running test: $$test"; \
		$(CXX) $(CFLAGS) $(EXTRA_CFLAGS) -c $(TEST_DIR)/$$test.cpp -o $(BUILD_DIR)/$$test.o; \
		$(CXX) $(BUILD_DIR)/$$test.o $(OBJS) -o $(BUILD_DIR)/$$test $(LDFLAGS); \
		if [ $$? -eq 0 ]; then \
			./$(BUILD_DIR)/$$test; \
		else \
			echo "Test $$test failed to compile."; \
		fi; \
	done

# Custom target to compile with maximum warnings
.PHONY: make_warning
make_warning: CFLAGS += $(EXTRA_CFLAGS)
make_warning: all

# Custom target to build simulation executable
simulation: $(OBJS) $(MAIN_OBJ)
	$(CXX) $(OBJS) $(MAIN_OBJ) -o $(BUILD_DIR)/simulation $(LDFLAGS)
