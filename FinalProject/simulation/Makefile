SRC_DIR     := src
INCLUDE_DIR := include
BUILD_DIR   := build
TEST_DIR    := test

CXX         := g++
CXXFLAGS    := -O2 -Wall -I$(INCLUDE_DIR)
EXTRA_FLAGS := -Wextra -pedantic -Wconversion -Wshadow

SRCS        := $(wildcard $(SRC_DIR)/*.cpp)
OBJS        := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

MAIN_SRC    := main.cpp
MAIN_OBJ    := $(BUILD_DIR)/main.o
EXECUTABLE  := $(BUILD_DIR)/main.exe

.PHONY: all
all: $(BUILD_DIR) $(OBJS)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compilation rule for main.cpp
$(MAIN_OBJ): $(MAIN_SRC) | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Linking rule for main.exe
$(EXECUTABLE): $(OBJS) $(MAIN_OBJ)
	@echo "Linking main.exe..."
	$(CXX) $(CXXFLAGS) $^ -o $@

.PHONY: run
run: $(EXECUTABLE)
	@echo "Executing main.exe..."
	./$(EXECUTABLE)

.PHONY: test
test:
	@if [ -z "$(TEST)" ]; then \
		echo "Error: No test specified. Usage: make test TEST=<test_name>"; \
		exit 1; \
	fi
	@echo "Compiling and running test: $(TEST)..."
	$(CXX) $(CXXFLAGS) $(EXTRA_FLAGS) $(TEST_DIR)/test_$(TEST).cpp $(OBJS) -o $(BUILD_DIR)/test_$(TEST)
	@echo "Executing test: $(TEST)..."
	./$(BUILD_DIR)/test_$(TEST)

.PHONY: clean
clean:
	@echo "Cleaning build directory..."
	rm -rf $(BUILD_DIR)