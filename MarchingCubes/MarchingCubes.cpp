// MarchingCubes.cpp : Defines the entry point for the application.

#include "MarchingCubes.h"
#include "mc_tables.h"

using namespace std;

#define ISOVALUE 0

int main(int argc, char** argv)
{
	// Read in the scalar field 
	std::ifstream dField(argv[1]);

	// Dimensions of the scalar field are first in the file
	long dimension;
	dField >> dimension;
	//DEBUG LINE
	//cout << "Dimension of field: " << dimension << endl;

	// Define the scalar field container size
	scalarField.resize(dimension, vector<vector<float>>(dimension, vector<float>(dimension)));
	float distance = 0;
	// Fill the scalar field container
	for (int x = 0; x < dimension; x++) {
		for (int y = 0; y < dimension; y++) {
			for (int z = 0; z < dimension; z++) {
				dField >> distance;
				scalarField[x][y][z] = distance;
			}
		}
	}
	dField.close();

	// Define vertex holding structures
	vector<vector<float>> vertices;
	int numVertices = -1;
	float vertexA[3];
	float vertexB[3];
	float distA;
	float distB;

	// March through the distance field
	for(int x = 0; x < dimension - 1; x++){
		for(int y = 0; y < dimension - 1; y++){
			for(int z = 0; z < dimension - 1; z++){
				// Get the state of the cube
				int caseNum = 0;
				for(int i = 0; i < 8; i++){
					// Check if inside the mesh
					if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] < ISOVALUE) {
						caseNum += pow(2, i);
					}
				}
				//DEBUG LINE
				//cout << "Vertex: "<< x << ", " << y << ", " << z << ", " << "Case:" << caseNum << endl;

				// Use lookup table to find case
				int numTriangles = triangleTable[caseNum][0];
				for(int i = 1; i < (numTriangles*3)+1; i++){
					// Edge to place vertex
					int edge = triangleTable[caseNum][i];
					//DEBUG LINE
					//cout << "Edge: " << triangleTable[caseNum][i] << endl;
					// Add a new vertex
					vertices.push_back(vector<float>());
					numVertices++;
					// Store the edge vertices
					for(int j = 0; j < 3; j++){
						vertexA[j] = vertPos[edgeTable[edge][0]][j];
						vertexB[j] = vertPos[edgeTable[edge][1]][j];
						//DEBUG LINE
						//vertices[numVertices].push_back(float(vertPos[edgeTable[edge][0]][j] + vertPos[edgeTable[edge][1]][j]) / 2.0);
					}
					// Add the current position to put into world co-ords
					vertexA[0] += x;
					vertexB[0] += x;
					vertexA[1] += y;
					vertexB[1] += y;
					vertexA[2] += z;
					vertexB[2] += z;
					// Get distance at the vertices from scalar fields
					distA = scalarField[vertexA[0]][vertexA[1]][vertexA[2]];
					distB = scalarField[vertexB[0]][vertexB[1]][vertexB[2]];
					// Get vertex position on the edge using linear interpolation
					for(int j = 0; j < 3; j++){
						// Store resulting triangle after doing linear interpolation
						vertices[numVertices].push_back(float(vertexA[j] + (ISOVALUE - distA) * (vertexB[j] - vertexA[j]) / (distB - distA)));
					}
					/*
					// DEBUG LINES
					cout << "Distance A: " << distA;
					cout << " Vertex A: " << vertexA[0] << ", " << vertexA[1] << ", " << vertexA[2] << endl;
					cout << "Distance B: " << distB;
					cout << " Vertex B: " << vertexB[0] << ", " << vertexB[1] << ", " << vertexB[2] << endl;
					cout << "Vertex P: " << vertices[numVertices][0] + x << ", " << vertices[numVertices][1] + y << ", " << vertices[numVertices][2] + z << endl;
					*/	
				}
			}
		}
	}

	// Write the output to triangle soup file
	ofstream out("out.tri");
	out <<  vertices.size() / 3 << endl;

	for(size_t i = 0; i < vertices.size(); i++){
		out << std::fixed << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << " " << endl;
	}

	out.close();

	return 0;
}