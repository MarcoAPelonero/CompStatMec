CXX = g++
CXXFLAGS = -I./include -std=c++11
LDFLAGS =

SRC = src/tensors.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = test

all: $(TARGET)

$(TARGET): $(OBJ) test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test.o: test.cpp
	$(CXX) $(CXXFLAGS) -c test.cpp

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) test.o $(TARGET)
