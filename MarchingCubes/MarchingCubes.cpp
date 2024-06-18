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
	cout << "Dimension of field: " << dimension << endl;

	// Define the scalar field container
	scalarField.resize(dimension, vector<vector<float>>(dimension, vector<float>(dimension)));
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

	// Define vertex holding structure
	vector<vector<float>> vertices;
	int numVertices = -1;

	// March through the distance field
	for(int x = 0; x < dimension - 1; x++){
		for(int y = 0; y < dimension - 1; y++){
			for(int z = 0; z < dimension - 1; z++){
				// Get the state of the cube
				int caseNum = 0;
				for(int i = 0; i < 8; i++){
					if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] < 0) {
						caseNum += pow(2, i);
					}
				}
				cout << "Vertex: "<< x << ", " << y << ", " << z << ", " << "Case:" << caseNum << endl;

				// Use lookup table to find case
				int numTriangles = triangleTable[caseNum][0];
				cout << "Number of triangles: " << numTriangles << endl;
				for(int i = 1; i < (numTriangles*3)+1; i++){
					// Edge to place vertex
					int edge = triangleTable[caseNum][i];
					cout << "Edge: " << triangleTable[caseNum][i] << " Vertex: ";
					vertices.push_back(vector<float>());
					numVertices++;
					// Vertex position in world co-ords
					for(int j = 0; j < 3; j++){
						// Store resulting triangle
						vertices[numVertices].push_back(float(vertPos[edgeTable[edge][0]][j] + vertPos[edgeTable[edge][1]][j]) / 2.0);
					}
					cout << vertices[numVertices][0] + x << ", " << vertices[numVertices][1] + y << ", " << vertices[numVertices][2] + z << endl;
				}
			}
		}
	}

	return 0;
}

int getCubeState(int x, int y, int z){
	int caseNum = 0;

	for(int i = 0; i < 8; i++){
		//cout << scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] << endl;
		if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] < 0) {
			caseNum += pow(2, i);
		}
	}

	cout << "Case:" << caseNum << endl;

	return caseNum;
}