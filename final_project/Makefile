CXXFLAGS = -I./include -std=c++11
LDFLAGS =

# Automatically find all .cpp files in the src directory
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)
TARGET = test

all: $(TARGET) clean-obj

$(TARGET): $(OBJ) test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

test.o: test.cpp
	$(CXX) $(CXXFLAGS) -MMD -c test.cpp

clean-obj:
	rm -f $(OBJ) test.o $(DEP)

clean:
	rm -f $(OBJ) test.o $(TARGET) $(DEP)

# Include the .d files to track header dependencies
-include $(DEP)