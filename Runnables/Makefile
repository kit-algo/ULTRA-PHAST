CC=g++
FLAGS=-std=c++17 -fopenmp -pipe
OPTIMIZATION=-march=native -O3
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wfatal-errors -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer
RELEASE=-ffast-math -ftree-vectorize -Wfatal-errors -DNDEBUG

all: BuildCH BuildCoreCH ComputeShortcuts RunCSAQueriesToBall RunCSAQueriesToStops RunCSAQueriesToVertices RunRAPTORQueriesToBall RunRAPTORQueriesToStops RunRAPTORQueriesToVertices

BuildCH:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o BuildCH BuildCH.cpp

BuildCoreCH:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o BuildCoreCH BuildCoreCH.cpp

ComputeShortcuts:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o ComputeShortcuts ComputeShortcuts.cpp

RunCSAQueriesToBall:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunCSAQueriesToBall RunCSAQueriesToBall.cpp
	
RunCSAQueriesToStops:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunCSAQueriesToStops RunCSAQueriesToStops.cpp

RunCSAQueriesToVertices:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunCSAQueriesToVertices RunCSAQueriesToVertices.cpp
	
RunRAPTORQueriesToBall:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunRAPTORQueriesToBall RunRAPTORQueriesToBall.cpp
	
RunRAPTORQueriesToStops:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunRAPTORQueriesToStops RunRAPTORQueriesToStops.cpp

RunRAPTORQueriesToVertices:
	$(CC) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o RunRAPTORQueriesToVertices RunRAPTORQueriesToVertices.cpp
