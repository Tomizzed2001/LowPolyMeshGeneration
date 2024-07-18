// MarchingCubes.cpp : Defines the entry point for the application.

#include "MarchingCubes.h"
#include "mc_tables.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define ISOVALUE 0

int main(int argc, char** argv){
	// Read in the scalar field 
	std::ifstream dField(argv[1]);

	// Dimensions of the scalar field are first in the file
	long xDimension, yDimension, zDimension;
	dField >> xDimension >> yDimension >> zDimension;

	// Define the scalar field container size
	scalarField.resize(xDimension, vector<vector<float>>(yDimension, vector<float>(zDimension)));
	float distance = 0;
	// Fill the scalar field container
	for (int x = 0; x < xDimension; x++) {
		for (int y = 0; y < yDimension; y++) {
			for (int z = 0; z < zDimension; z++) {
				dField >> distance;
				scalarField[x][y][z] = distance;
			}
		}
	}
	dField.close();

	// Define vertex holding structures
	vector<glm::vec3> triSoup;
	glm::vec3 vertexA;
	glm::vec3 vertexB;
	float distA;
	float distB;

	//xDimension = 20;
	//yDimension = 14;
	//zDimension = 14;


	// March through the distance field
	// TEST CASE IS 16 15 18
	for(int x = 0; x < xDimension - 1; x++){
		for(int y = 0; y < yDimension - 1; y++){
			for(int z = 0; z < zDimension - 1; z++){
				//cout << "New Cube" << endl;
				// Get the state of the cube
				int caseNum = 0;
				for(int i = 0; i < 8; i++){
					// Check if inside the mesh
					if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] >= ISOVALUE) {
						caseNum += pow(2, i);
					}
					//cout << "Point: " << x+vertPos[i][0] << "," << y+vertPos[i][1] << "," << z+vertPos[i][2] << " Distance: " << scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] << " CaseNum: " << caseNum <<  endl;
				}
				//cout << "Case Num: " << caseNum << endl;

				// Use lookup table to find case
				int numTriangles = triangleTable[caseNum][0];
				//cout << "Number of triangles: " << numTriangles << endl;
				for(int i = 1; i < (numTriangles*3)+1; i++){
					// Edge to place vertex
					int edge = triangleTable[caseNum][i];
					//cout << "Edge: " << edge << endl;

					// Store the edge vertices
					for(int j = 0; j < 3; j++){
						vertexA[j] = vertPos[edgeTable[edge][0]][j];
						vertexB[j] = vertPos[edgeTable[edge][1]][j];
					}
					//cout << "VertexA: " << glm::to_string(vertexA) << endl;
					//cout << "VertexB: " << glm::to_string(vertexB) << endl;

					// Add the current position to put into world co-ords
					vertexA = vertexA + glm::vec3(x,y,z);
					vertexB = vertexB + glm::vec3(x,y,z);

					// Get distance at the vertices from scalar fields
					distA = scalarField[vertexA[0]][vertexA[1]][vertexA[2]];
					distB = scalarField[vertexB[0]][vertexB[1]][vertexB[2]];

					// Get vertex position on the edge using linear interpolation
					//triSoup.push_back(vertexA + (ISOVALUE - distA) * (vertexB - vertexA) / (distB - distA));
					triSoup.push_back((vertexA + vertexB)/float(2.0));
				}
				//break; 
			}
			//break;
		}
		//break;
	} 

	/* Output to triangle soup .tri file
	// Write the output to triangle soup file
	ofstream outOld("out.tri");
	outOld <<  triSoup.size() / 3 << endl;

	for(size_t i = 0; i < triSoup.size(); i++){
		outOld << std::fixed << triSoup[i][0] << " " << triSoup[i][1] << " " << triSoup[i][2] << " " << endl;
	}

	outOld.close();
	*/

	// Build the directed edge data structure from the triangle soup.
	for(size_t vID = 0; vID < triSoup.size(); vID++){
		// Get the vertex index and place if not found.
		int vertex = getVertexID(triSoup[vID]);
		// Add the vertex id to the face
		faces.push_back(vertex);
	}

	// Index storage variables
	int v0, v1;

	otherHalf.resize(faces.size(),-1);

	// Get the other half
	for(size_t eID = 0; eID < faces.size(); eID++){
		// Since faces form the edges check if it needs to search for 2 to 0 edge
		if (eID % 3 == 2){
			v0 = faces[eID];
			v1 = faces[eID-2];
		}
		else{
			v0 = faces[eID];
			v1 = faces[eID+1];
		}

		// Look for v1 to v0
		for(size_t fID = 0; fID < faces.size(); fID+=3){
			if(faces[fID] == v1 && faces[fID+1] == v0){
				otherHalf[eID] = fID + 0;
				break;
			}
			else if(faces[fID+1] == v1 && faces[fID+2] == v0){
				otherHalf[eID] = fID + 1;
				break;
			}
			else if(faces[fID+2] == v1 && faces[fID] == v0){
				otherHalf[eID] = fID + 2;
				break;
			}
		}
	}

	// Write the output to a directed edge file format
	ofstream out("out.diredge");
	for(size_t i = 0; i < vertices.size(); i++){
		out << "v " << std::fixed << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl;
	}
	for(size_t i = 0; i < firstDirectedEdges.size(); i++){
		out << "fde " << firstDirectedEdges[i] << endl;
	} 
	for(size_t i = 0; i < faces.size(); i+=3){
		out << "f " << faces[i] << " " << faces[i+1] << " " << faces[i+2] << endl;
	}
	for(size_t i = 0; i < otherHalf.size(); i++){
		out << "oh " << otherHalf[i] << endl;
	}

	out.close();

	return 0;
}

int getVertexID(glm::vec3 vertex){
	// Look to see if the vertex already exists
	for(size_t i = 0; i < vertices.size(); i++){
		if(vertex == vertices[i]){
			return i;
		}
	}
	// Place the vertex in the vertex array and return the index
	vertices.push_back(vertex);
	firstDirectedEdges.push_back(faces.size());
	return vertices.size() - 1;
}