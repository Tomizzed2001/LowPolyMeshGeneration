CC = g++
OBJ = DistanceFields.o
CFLAGS = -c -Wall -Iinclude

DistanceFields: $(OBJ)
	$(CC)  $(OBJ) -o DistanceFields

DistanceFields.o: DistanceFields.cpp DistanceFields.h
	$(CC) $(CFLAGS) DistanceFields.cpp

clean:
	rm *o DistanceFields
