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

    vector<glm::vec3> vertices;
    vector<glm::vec3> normals;
    vector<glm::vec3> faces;
    float minPoint[3] = {0}, maxPoint[3] = {0};

	// Read in the .obj file
	ifstream dField(argv[1]);
    if (dField.is_open()){
    
        // Initialise helper variables
        string newLine;
        string ss;
        string vtn;

        // Loop through .obj file line by line
        while (getline(dField, newLine)) {
            stringstream line(newLine);
            line >> ss;
            // If vertex normal
            if(ss == "vn"){
                glm::vec3 n;
                line >> n.x >> n.y >> n.z;
                normals.push_back(n);
            }
            // If vertex
            else if(ss == "v"){
                glm::vec3 p;
                // Format is x y z as a point
                for(int i = 0; i < 3; i++){
                    line >> p[i];
                    // Check if its a minimum or maximum
                    if(p[i] < minPoint[i]){
                        minPoint[i] = p[i];
                    }
                    else if(p[i] > maxPoint[i]){
                        maxPoint[i] = p[i];
                    }
                }
                vertices.push_back(p);
            }
            
            // If face
            else if(ss == "f"){
                glm::vec3 f;
                // Format is v/vt/vn so split these into the face vector as v1,n1,v2,n2,v3,n3
                for(int i = 0; i < 3; i++){
                    line >> vtn;
                    stringstream ssFace(vtn);
                    // Get Vertex
                    getline(ssFace, vtn, '/');
                    f[i] = stoi(vtn) - 1;
                    // Ignore the texture coords
                    getline(ssFace, vtn, '/');
                    // Get Normal
                    getline(ssFace, vtn, '/');
                }
                faces.push_back(f);
                //cout << "Face: " << glm::to_string(faces[11]) << endl;
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
        cout << glm::to_string(faces[faceID]) << endl;
        transforms.push_back(getTransformMatrix(vertices[faces[faceID].x], vertices[faces[faceID].y], vertices[faces[faceID].z]));
    }

    // Loop through all points in the scalar field and get the distance 
    // from them to the closest triangle
    for (float xID = 0, x = minPoint[0]; xID < xSize; xID++, x+=GRID_SIZE) {
		for (float yID = 0, y = minPoint[1]; yID < ySize; yID++, y+=GRID_SIZE) {
			for (float zID = 0, z = minPoint[2]; zID < zSize; zID++, z+=GRID_SIZE) {
                
            }
        }
    }

    // Testing stuff (It seems that it flips the distance so - is on the outside? Requires further testing to be
    // sure but test with real distance and can do the paper calculations too).
    glm::vec3 A = vertices[1];
    glm::vec3 B = vertices[7];
    glm::vec3 C = vertices[5];
    glm::vec3 N = glm::cross(glm::normalize(B-A), glm::normalize(C-A));

    glm::mat4 R = getTransformMatrix(A, B, C);

    A = R * glm::vec4(A,1);
    B = R * glm::vec4(B,1);
    C = R * glm::vec4(C,1);
    N = R * glm::vec4(N,1);

    cout << "A: " << A.x << ", " << A.y << ", " << A.z << endl;
    cout << "B: " << B.x << ", " << B.y << ", " << B.z << endl;
    cout << "C: " << C.x << ", " << C.y << ", " << C.z << endl;
    cout << "N: " << N.x << ", " << N.y << ", " << N.z << endl;
    
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