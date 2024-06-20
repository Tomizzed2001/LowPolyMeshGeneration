CC = g++
OBJ = DistanceFields.o
CFLAGS = -c -Wall -Iinclude

MarchingCubes: $(OBJ)
	$(CC)  $(OBJ) -o DistanceFields

MarchingCubes.o: DistanceFields.cpp DistanceFields.h
	$(CC) $(CFLAGS) DistanceFields.cpp

clean:
	rm *o DistanceFields
