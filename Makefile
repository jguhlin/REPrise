CXX := g++
CXXFLAGS := -std=c++11 -fopenmp -Ofast -Wall -Wextra -pedantic

all: REPrise

REPrise: REPrise.cpp sais_long.c
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm REPrise
