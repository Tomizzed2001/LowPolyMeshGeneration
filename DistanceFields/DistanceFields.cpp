// DistanceFields.cpp : Defines the entry point for the application.

#include "DistanceFields.h"

using namespace std;

int main(int argc, char** argv)
{
    vector<vector<float>> vertices;
    int numVertices = -1;
    vector<vector<float>> normals;
    int numNormals = -1;
    vector<vector<int>> faces;
    int numFaces = -1;

	// Read in the .obj file
	ifstream dField(argv[1]);
    if (dField.is_open()){
    
        // Initialise helper variables
        string newLine;
        string ss;
        string vtn;
        float point;
        // Loop through .obj file line by line
        while (getline(dField, newLine)) {
            stringstream line(newLine);
            line >> ss;
            // If vertex normal
            if(ss == "vn"){
                normals.push_back(vector<float>());
                numNormals++;
                // Format is x y z as a direction
                for(int i = 0; i < 3; i++){
                    line >> point;
                    normals[numNormals].push_back(point);
                }
                cout << "Normal: " << normals[numNormals][0] << ", " << normals[numNormals][1] << ", " << normals[numNormals][2] << endl;
            }
            // If vertex
            else if(ss == "v"){
                vertices.push_back(vector<float>());
                numVertices++;
                // Format is x y z as a point
                for(int i = 0; i < 3; i++){
                    line >> point;
                    vertices[numVertices].push_back(point);
                }
                cout << "Vertex: " << vertices[numVertices][0] << ", " << vertices[numVertices][1] << ", " << vertices[numVertices][2] << endl;
            }
            // If face
            else if(ss == "f"){
                faces.push_back(vector<int>());
                numFaces++;
                // Format is v/vt/vn so split these into the face vector as v1,n1,v2,n2,v3,n3
                for(int i = 0; i < 3; i++){
                    line >> vtn;
                    stringstream ssFace(vtn);
                    // Get Vertex
                    getline(ssFace, vtn, '/');
                    faces[numFaces].push_back(stoi(vtn));
                    // Ignore the texture coords
                    getline(ssFace, vtn, '/');
                    // Get Normal
                    getline(ssFace, vtn, '/');
                    faces[numFaces].push_back(stoi(vtn));
                }
                cout << "Face: " << faces[numFaces][0] << ", " << faces[numFaces][2] << ", " << faces[numFaces][4] << endl;
            }
        }
    }
    dField.close();

    // Set up the scalar field
    
	return 0;
}