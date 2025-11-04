CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET = fish_simulation
SOURCES = main.cpp

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) -lm

clean:
	rm -f $(TARGET)

.PHONY: clean