// DistanceFields.cpp : Defines the entry point for the application.

#include "DistanceFields.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define GRID_SIZE 0.5

int main(int argc, char** argv)
{
    //Check the file is a .obj file
    std::string s = argv[1];
    if (!s.find(".obj")){
        std::cout << "Program only takes input of a .obj file" << std::endl;
        return 0;
    }

    vector<vector<float>> vertices;
    int numVertices = -1;
    vector<vector<float>> normals;
    int numNormals = -1;
    vector<vector<int>> faces;
    int numFaces = -1;
    float minPoint[3] = {0}, maxPoint[3] = {0};

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
                //cout << "Normal: " << normals[numNormals][0] << ", " << normals[numNormals][1] << ", " << normals[numNormals][2] << endl;
            }
            // If vertex
            else if(ss == "v"){
                vertices.push_back(vector<float>());
                numVertices++;
                // Format is x y z as a point
                for(int i = 0; i < 3; i++){
                    line >> point;
                    vertices[numVertices].push_back(point);
                    // Check if its a minimum or maximum
                    if(point < minPoint[i]){
                        minPoint[i] = point;
                    }
                    else if(point > maxPoint[i]){
                        maxPoint[i] = point;
                    }
                }
                //cout << "Vertex: " << vertices[numVertices][0] << ", " << vertices[numVertices][1] << ", " << vertices[numVertices][2] << endl;
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
                //cout << "Face: " << faces[numFaces][0] << ", " << faces[numFaces][2] << ", " << faces[numFaces][4] << endl;
            }
        }
    }
    dField.close();

    
    //Round the points so the grid fits nicely
    cout << "Min Point: " << minPoint[0] << ", " << minPoint[1] << ", " << minPoint[2] << endl;
    cout << "Max Point: " << maxPoint[0] << ", " << maxPoint[1] << ", " << maxPoint[2] << endl;
    fitToGrid(minPoint, true);
    fitToGrid(maxPoint, false);
    cout << "Min Point: " << minPoint[0] << ", " << minPoint[1] << ", " << minPoint[2] << endl;
    cout << "Max Point: " << maxPoint[0] << ", " << maxPoint[1] << ", " << maxPoint[2] << endl;

    // Set up the scalar field
    int xSize = ((maxPoint[0] - minPoint[0]) / GRID_SIZE) + 1;
    int ySize = ((maxPoint[1] - minPoint[1]) / GRID_SIZE) + 1;
    int zSize = ((maxPoint[2] - minPoint[2]) / GRID_SIZE) + 1;
    std::vector<std::vector<std::vector<float>>> scalarField(xSize, vector<vector<float>>(ySize, vector<float>(zSize)));

    // Get the transform matrices for each triangle
    vector<glm::mat4> transforms;
    for(size_t faceID = 0; faceID < faces.size(); faceID++){

    }



    // Loop through all points in the scalar field and get the distance 
    // from them to the closest triangle
    for (float xID = 0, x = minPoint[0]; xID < xSize; xID++, x+=GRID_SIZE) {
		for (float yID = 0, y = minPoint[1]; yID < ySize; yID++, y+=GRID_SIZE) {
			for (float zID = 0, z = minPoint[2]; zID < zSize; zID++, z+=GRID_SIZE) {
                
            }
        }
    }

    glm::vec3 A = glm::vec3(3,1,-4);
    glm::vec3 B = glm::vec3(6,1,0);
    glm::vec3 C = glm::vec3(4,5,6);

    glm::mat4 R = getTransformMatrix(A, B, C);

    A = R * glm::vec4(A,1);
    B = R * glm::vec4(B,1);
    C = R * glm::vec4(C,1);

    cout << "A: " << A.x << ", " << A.y << ", " << A.z << endl;
    cout << "B: " << B.x << ", " << B.y << ", " << B.z << endl;
    cout << "C: " << C.x << ", " << C.y << ", " << C.z << endl;
    
	return 0;
}

void fitToGrid(float *point, bool isMin){
    if(isMin){
        point[0] = -(abs(point[0]) + (GRID_SIZE - fmod(abs(point[0]), GRID_SIZE))) - GRID_SIZE;
        point[1] = -(abs(point[1]) + (GRID_SIZE - fmod(abs(point[1]), GRID_SIZE))) - GRID_SIZE;
        point[2] = -(abs(point[2]) + (GRID_SIZE - fmod(abs(point[2]), GRID_SIZE))) - GRID_SIZE;  
        return;
    }
    else{
        point[0] = point[0] + (GRID_SIZE - fmod(point[0], GRID_SIZE)) + GRID_SIZE;
        point[1] = point[1] + (GRID_SIZE - fmod(point[1], GRID_SIZE)) + GRID_SIZE;
        point[2] = point[2] + (GRID_SIZE - fmod(point[2], GRID_SIZE)) + GRID_SIZE;  
        return;
    }
}

glm::mat4 getTransformMatrix(glm::vec3 A, glm::vec3 B, glm::vec3 C){
    glm::vec4 translate = glm::vec4(-A, 1);

    glm::vec3 U = glm::normalize(B - A);
    glm::vec3 W = glm::normalize(glm::cross(U, C - A));
    glm::vec3 V = glm::cross(U,W);

    glm::mat3 RotMatrix = glm::transpose(glm::mat3(W,V,U));

    glm::mat4 R = glm::mat4(RotMatrix);
    R[3][3] = 1;

    glm::mat4 T = glm::mat4(1.f);
    T[3] = translate;

    return R * T;

}