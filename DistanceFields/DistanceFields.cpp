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

    vector<glm::vec3> vertices; // x y z
    vector<glm::vec3> normals;  // x y z directions
    vector<glm::vec3> faces;    // vert vert vert
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
    transformedFaces.resize(faces.size()+1, vector<glm::vec2>(6));
    for(size_t faceID = 0; faceID < faces.size(); faceID++){
        transforms.push_back(getTransformMatrix(vertices[faces[faceID].x], vertices[faces[faceID].y], vertices[faces[faceID].z]));
        preComputeFace(faceID, transforms.back(), vertices[faces[faceID].x], vertices[faces[faceID].y], vertices[faces[faceID].z]);
    }

    // Loop through all points in the scalar field and get the distance 
    // from them to the closest triangle
    for (float xID = 0, x = minPoint[0]; xID < xSize; xID++, x+=GRID_SIZE) {
		for (float yID = 0, y = minPoint[1]; yID < ySize; yID++, y+=GRID_SIZE) {
			for (float zID = 0, z = minPoint[2]; zID < zSize; zID++, z+=GRID_SIZE) {
                float closestDistance = 1000;
                for(size_t faceID = 0; faceID < faces.size(); faceID++){
                    /*
                    float distance = distToTriangle(
                        transform,
                        glm::vec3(x,y,z),
                        vertices[faces[faceID].x],
                        vertices[faces[faceID].y], 
                        vertices[faces[faceID].z]
                    );
                    if(abs(distance) < closestDistance){
                        closestDistance = distance;
                    }
                    */
                }
                scalarField[xID][yID][zID] = closestDistance;
            }
        }
    }

    // Dump the scalar field to a .txt file

    // Testing stuff (It seems that it flips the distance so - is on the outside? Requires further testing to be
    // sure but test with real distance and can do the paper calculations too).
    glm::mat4 T = glm::mat4(1);
    transforms[0] = T;
    glm::vec3 A = glm::vec3(0,0,0);
    glm::vec3 B = glm::vec3(0,0,3);
    glm::vec3 C = glm::vec3(0,3,0);
    glm::vec3 P = glm::vec3(5,1,1);

    preComputeFace(0, T, A, B, C);
    cout << distToTriangle(0, P) << endl;
    
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

void preComputeFace(int faceID, glm::mat4 transform, glm::vec3 A, glm::vec3 B, glm::vec3 C){
    // Transform the original points for 0 - 2
    glm::vec4 tA = transform * glm::vec4(A,1);
    transformedFaces[faceID][0] = glm::vec2(tA.z, tA.y);
    glm::vec4 tB = transform * glm::vec4(B,1);
    transformedFaces[faceID][1] = glm::vec2(tB.z, tB.y);
    glm::vec4 tC = transform * glm::vec4(C,1);
    transformedFaces[faceID][2] = glm::vec2(tC.z, tC.y);

    // Get the normals
    glm::vec2 AB = transformedFaces[faceID][1] - transformedFaces[faceID][0];
    transformedFaces[faceID][3] = glm::vec2(AB.y, -AB.x);
    glm::vec2 BC = transformedFaces[faceID][2] - transformedFaces[faceID][1];
    transformedFaces[faceID][4] = glm::vec2(BC.y, -BC.x);
    glm::vec2 CA = transformedFaces[faceID][0] - transformedFaces[faceID][2];
    transformedFaces[faceID][5] = glm::vec2(CA.y, -CA.x);
}

float distToTriangle(int faceID, glm::vec3 P){
    // Transform P into 2D plane of triangle
    P = transforms[faceID] * glm::vec4(P,1);
    glm::vec2 tP = glm::vec2(P.z, P.y);
    glm::vec2 tA = transformedFaces[faceID][0];
    glm::vec2 tB = transformedFaces[faceID][1];
    glm::vec2 tC = transformedFaces[faceID][2];

    // Determine if point is within the bounds of the triangle using the half plane test
    float AB = (tB.y - tA.y) * tP.x + (tA.x - tB.x) * tP.y + (tB.x * tA.y - tA.x * tB.y);
    float BC = (tC.y - tB.y) * tP.x + (tB.x - tC.x) * tP.y + (tC.x * tB.y - tB.x * tC.y);
    float CA = (tA.y - tC.y) * tP.x + (tC.x - tA.x) * tP.y + (tA.x * tC.y - tC.x * tA.y);
    if(AB <= 0 && BC <= 0 && CA <= 0){
        return P.x;
    }

    // Point must be outside the triangle so check to see if there is a closest vertex
    // Add the normals onto the required points and do more 2D tests oyeah


    return 0;
    
}