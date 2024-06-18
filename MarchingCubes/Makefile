CC = g++
OBJ = MarchingCubes.o
CFLAGS = -c -Wall -Iinclude

MarchingCubes: $(OBJ)
	$(CC)  $(OBJ) -o MarchingCubes

MarchingCubes.o: MarchingCubes.cpp MarchingCubes.h
	$(CC) $(CFLAGS) MarchingCubes.cpp

clean:
	rm *o marchingCubes
