// MarchingCubes.cpp : Defines the entry point for the application.
//

#include "MarchingCubes.h"
#include "mc_tables.h"

using namespace std;

int main(int argc, char** argv)
{
	// Read in the scalar field 
	std::ifstream dField(argv[1]);

	// Dimensions of the scalar field are first in the file
	long dimension;
	dField >> dimension;
	cout << dimension << endl;

	// Define the scalar field container
	std::vector<std::vector<std::vector<float>>> scalarField(dimension, vector<vector<float>>(dimension, vector<float>(dimension)));
	int distance = 0;
	for (int x = 0; x < dimension; x++) {
		for (int y = 0; y < dimension; y++) {
			for (int z = 0; z < dimension; z++) {
				dField >> distance;
				scalarField[x][y][z] = distance;
			}
		}
	}

	dField.close();

	// March through the distance field
    cout << edgeTable[2][1] << endl;

	// Get the state of the cube

	// Use lookup table

	// Store resulting triangles


	return 0;
}
