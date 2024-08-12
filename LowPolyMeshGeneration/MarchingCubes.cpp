// MarchingCubes.cpp : Defines the entry point for the application.

#include "MarchingCubes.h"
#include "mc_tables_correct.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define ISOVALUE 0.00005

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

	// Get the optional inputs from the user
    bool outputOBJ = false;
    bool answered = false;
    string answer;
    while (!answered)
    {
        cout << "Would you like to output an OBJ file? (y/n)";
        cin >> answer;
        if(answer == "y" || answer == "Y"){
			outputOBJ = true;
            answered = true;
        }
        if(answer == "n" || answer == "N"){
            answered = true;
        }
    }

	// TIMER START
    auto start1 = std::chrono::high_resolution_clock::now();

	array<int, 3> empty = {-1,-1,-1};
	cubeEdges.resize(xDimension-1, vector<vector<std::array<int,3>>>(yDimension-1, vector<std::array<int,3>>(zDimension-1, empty)));

	// Define vertex holding structures
	glm::vec3 vertexA;
	glm::vec3 vertexB;
	float distA;
	float distB;

	cout << "Beginning Marching Cubes" << endl;
	// March through the distance field
	for(int x = 0; x < xDimension - 1; x++){
		for(int y = 0; y < yDimension - 1; y++){
			for(int z = 0; z < zDimension - 1; z++){
				// Get the state of the cube
				int caseNum = 0;
				for(int i = 0; i < 8; i++){
					// Check if inside the mesh
					if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] == 1000){
						caseNum = 0;
						break;
					}
					if (scalarField[x+vertPos[i][0]][y+vertPos[i][1]][z+vertPos[i][2]] >= ISOVALUE) {
						caseNum += pow(2, i);
					}
				}

				// Use lookup table to find case
				int numTriangles = triangleTable[caseNum][0];
				for(int i = 1; i < (numTriangles*3)+1; i++){
					// Edge to place vertex
					int edge = triangleTable[caseNum][i];

					// Store the edge vertices
					for(int j = 0; j < 3; j++){
						vertexA[j] = vertPos[edgeTable[edge][0]][j];
						vertexB[j] = vertPos[edgeTable[edge][1]][j];
					}

					// Add the current position to put into world co-ords
					vertexA = vertexA + glm::vec3(x,y,z);
					vertexB = vertexB + glm::vec3(x,y,z);

					// Get distance at the vertices from scalar fields
					distA = scalarField[vertexA[0]][vertexA[1]][vertexA[2]];
					distB = scalarField[vertexB[0]][vertexB[1]][vertexB[2]];

					// Get the edge of the cube in the cube edge table to see if the edge already has a vertex
					int currentEdgeValue = cubeEdges[cubeEdgeLookUp[edge][0] + x][cubeEdgeLookUp[edge][1] + y][cubeEdgeLookUp[edge][2] + z][cubeEdgeLookUp[edge][3]];

					if(currentEdgeValue == -1){
						// Nothing placed yet so place it
						triSoup.push_back(vertexA + float(ISOVALUE - distA) * (vertexB - vertexA) / (distB - distA));
						cubeEdges[cubeEdgeLookUp[edge][0] + x][cubeEdgeLookUp[edge][1] + y][cubeEdgeLookUp[edge][2] + z][cubeEdgeLookUp[edge][3]] = triSoup.size()-1;
					}
					else{
						triSoup.push_back(triSoup[currentEdgeValue]);
					}
				}
			}
		}
	} 

	cout << "Generating data structure" << endl;
	// Build the directed edge data structure from the triangle soup.
	for(size_t vID = 0; vID < triSoup.size(); vID++){
		// Get the vertex index and place if not found.
		int vertex = getVertexID(triSoup[vID]);
		// Add the vertex id to the face
		mesh.edges.push_back(vertex);
	}

	// Index storage variables
	int v0, v1;

	mesh.otherhalves.resize(mesh.edges.size(),-1);

	// Get the other half
	for(size_t eID = 0; eID < mesh.edges.size(); eID++){
		// Since faces form the edges check if it needs to search for 2 to 0 edge
		if (eID % 3 == 2){
			v0 = mesh.edges[eID];
			v1 = mesh.edges[eID-2];
		}
		else{
			v0 = mesh.edges[eID];
			v1 = mesh.edges[eID+1];
		}

		// Look for v1 to v0
		for(size_t fID = 0; fID < mesh.edges.size(); fID+=3){
			if(mesh.edges[fID] == v1 && mesh.edges[fID+1] == v0){
				mesh.otherhalves[eID] = fID + 0;
				break;
			}
			else if(mesh.edges[fID+1] == v1 && mesh.edges[fID+2] == v0){
				mesh.otherhalves[eID] = fID + 1;
				break;
			}
			else if(mesh.edges[fID+2] == v1 && mesh.edges[fID] == v0){
				mesh.otherhalves[eID] = fID + 2;
				break;
			}
		}
	}

	// TIMER 1 END
    auto end1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    cout << "Time taken: " << float(duration1.count()) / 1000000.0 << " seconds." << endl;

	// Write the output to a directed edge file format
	ofstream out("isosurface.diredge");
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
	for(size_t i = 0; i < mesh.edges.size(); i+=3){
		out << "f " << mesh.edges[i] << " " << mesh.edges[i+1] << " " << mesh.edges[i+2] << endl;
	}
	for(size_t i = 0; i < mesh.otherhalves.size(); i++){
		out << "oh " << mesh.otherhalves[i] << endl;
	}

	out.close();

	if(outputOBJ){
		outputToObject();
	}

	return 0;
}

int getVertexID(glm::vec3 vertex){
	// Look to see if the vertex already exists
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		if(vertex == mesh.vertices[i]){
			return i;
		}
	}
	// Place the vertex in the vertex array and return the index
	mesh.vertices.push_back(vertex);
	//mesh.firstDirectedEdges.push_back(faces.size());
	return mesh.vertices.size() - 1;
}

void outputToObject(){
    std::string var = "isosurface.obj";
    ofstream out(var);
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
	for(size_t i = 0; i < mesh.edges.size(); i+=3){
		out << "f " 
        << mesh.edges[i] + 1 << "//" << " "
        << mesh.edges[i+1] + 1 << "//" << " "
        << mesh.edges[i+2] + 1 << "//"
        << endl;
	}
	out.close();
}