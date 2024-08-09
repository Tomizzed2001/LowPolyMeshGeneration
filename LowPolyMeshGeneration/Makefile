CC = g++
CFLAGS = -c -Wall -Iinclude -fopenmp

all: DistanceFields MarchingCubes Simplification ManifoldChecker

DistanceFields: DistanceFields.o
	$(CC)  -fopenmp DistanceFields.o -o DistanceFields

DistanceFields.o: DistanceFields.cpp DistanceFields.h
	$(CC) $(CFLAGS) DistanceFields.cpp

MarchingCubes: MarchingCubes.o
	$(CC)  MarchingCubes.o -o MarchingCubes

MarchingCubes.o: MarchingCubes.cpp MarchingCubes.h
	$(CC) $(CFLAGS) MarchingCubes.cpp

Simplification: Simplification.o
	$(CC)  -fopenmp Simplification.o -o Simplification

Simplification.o: Simplification.cpp Simplification.h
	$(CC) $(CFLAGS) Simplification.cpp

ManifoldChecker: ManifoldChecker.o
	$(CC)  ManifoldChecker.o -o ManifoldChecker

ManifoldChecker.o: ManifoldChecker.cpp ManifoldChecker.h
	$(CC) $(CFLAGS) ManifoldChecker.cpp
