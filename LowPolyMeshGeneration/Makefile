CC = g++
CFLAGS = -c -Wall -Iinclude

all: DistanceFields MarchingCubes

DistanceFields: DistanceFields.o
	$(CC)  DistanceFields.o -o DistanceFields

DistanceFields.o: DistanceFields.cpp DistanceFields.h
	$(CC) $(CFLAGS) DistanceFields.cpp

MarchingCubes: MarchingCubes.o
	$(CC)  MarchingCubes.o -o MarchingCubes

MarchingCubes.o: MarchingCubes.cpp MarchingCubes.h
	$(CC) $(CFLAGS) MarchingCubes.cpp
