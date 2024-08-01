// DistanceFields.cpp : Defines the entry point for the application.

#include "DistanceFields.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define GRID_SIZE_MAX 64.0

int main(int argc, char** argv)
{
    //Check the file is a .obj file
    std::string s = argv[1];
    if (!s.find(".obj")){
        std::cout << "Program only takes input of a .obj file" << std::endl;
        return 0;
    }

    vector<glm::vec3> normals;  // x y z directions
    float minPoint[3] = {0}, maxPoint[3] = {0};
    bool firstVertex = true;

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
            // If vertex
            if(ss == "v"){
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
                inputMesh.vertices.push_back(p);
                if(firstVertex){
                    minPoint[0] = p.x;
                    maxPoint[0] = p.x;
                    minPoint[1] = p.y;
                    maxPoint[1] = p.y;
                    minPoint[2] = p.z;
                    maxPoint[2] = p.z;
                    firstVertex = false;
                }
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
                inputMesh.faces.push_back(f);
                //cout << "Face: " << glm::to_string(faces[11]) << endl;
            }
        }
    }
    dField.close();
    
    // TIMER START
    auto start1 = std::chrono::high_resolution_clock::now();

    cout << "Calculating face normals" << endl;
    // Calculate face normals from scratch
    fNormals.resize(inputMesh.faces.size());
    #pragma omp parallel for
    for(size_t fID = 0; fID < inputMesh.faces.size(); fID++){
        glm::vec3 AB = inputMesh.vertices[inputMesh.faces[fID].y] - inputMesh.vertices[inputMesh.faces[fID].x];
        glm::vec3 AC = inputMesh.vertices[inputMesh.faces[fID].z] - inputMesh.vertices[inputMesh.faces[fID].x];
        fNormals[fID] = glm::normalize(glm::cross(AB, AC));
    }

    cout << "Calculating vertex normals" << endl;
    // Calculate the vertex normals from scratch
    inputMesh.vNormals.resize(inputMesh.vertices.size());
    #pragma omp parallel for
    for(size_t vID = 0; vID < inputMesh.vertices.size(); vID++){
        inputMesh.vNormals[vID] = getVertexNormal(vID);
    }

    cout << "Calculating edge normals" << endl;
    // Calculate the edge normals from scratch
    eNormals.resize(inputMesh.faces.size()*3);
    #pragma omp parallel for
    for(size_t fID = 0; fID < inputMesh.faces.size(); fID++){
        // AB edge
        eNormals[fID * 3] = glm::normalize((findOtherHalfNormal(inputMesh.faces[fID].x, inputMesh.faces[fID].y) + fNormals[fID]) / float(2));
        // BC edge
        eNormals[fID * 3 + 1] = glm::normalize((findOtherHalfNormal(inputMesh.faces[fID].y, inputMesh.faces[fID].z) + fNormals[fID]) / float(2));
        // CA edge
        eNormals[fID * 3 + 2] = glm::normalize((findOtherHalfNormal(inputMesh.faces[fID].z, inputMesh.faces[fID].x) + fNormals[fID]) / float(2));
    }

    sizeOfGrid = max((maxPoint[0] - minPoint[0]), max((maxPoint[1] - minPoint[1]),(maxPoint[2] - minPoint[2]))) / GRID_SIZE_MAX;

    cout << sizeOfGrid << endl;

    //Round the points so the grid fits nicely
    fitToGrid(minPoint, true);
    fitToGrid(maxPoint, false);

    // Set up the scalar field
    int xSize = round(((maxPoint[0] - minPoint[0]) / sizeOfGrid) + 1);
    int ySize = round(((maxPoint[1] - minPoint[1]) / sizeOfGrid) + 1);
    int zSize = round(((maxPoint[2] - minPoint[2]) / sizeOfGrid) + 1);
    std::vector<std::vector<std::vector<float>>> scalarField(xSize, vector<vector<float>>(ySize, vector<float>(zSize)));


    cout << "Beginning distance calculations" << endl;
    
    // Loop through all points in the scalar field and get the distance 
    // from them to the closest triangle
    float xStart = minPoint[0];
    #pragma omp parallel for
    for (int xID = 0; xID < xSize; xID++) {
        float x = xStart + (xID * sizeOfGrid);
        float yStart = minPoint[1];
		for (int yID = 0; yID < ySize; yID++) {
            float y = yStart + (yID * sizeOfGrid);
            float zStart = minPoint[2];
			for (int zID = 0; zID < zSize; zID++) {
                float z = zStart + (zID * sizeOfGrid);
                float closestDistance = 1000;
                for(size_t faceID = 0; faceID < inputMesh.faces.size(); faceID++){
                    float distance = distanceToTriangle(
                        glm::vec3(x,y,z),               // Point
                        inputMesh.vertices[inputMesh.faces[faceID].x],      // Vertex A
                        inputMesh.vertices[inputMesh.faces[faceID].y],      // Vertex B
                        inputMesh.vertices[inputMesh.faces[faceID].z],      // Vertex C
                        faceID                          // Face ID
                    );

                    if(abs(distance) < abs(closestDistance)){
                        closestDistance = distance;
                    }
                }
                // Remove any "negative" zeros
                if(closestDistance <= 0.0001 && closestDistance >= -0.0001){
                    closestDistance = 0;
                }
                scalarField[xID][yID][zID] = closestDistance;
            }
        }
    }

    // TIMER 1 END
    auto end1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    cout << "Time taken: " << float(duration1.count()) / 1000000.0 << " seconds." << endl;

    // Dump the scalar field to a .txt file
    std::ofstream out("scalarField.txt");
    out << xSize << " " << ySize << " " << zSize << endl;
    for (int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			for (int z = 0; z < zSize; z++) {
                if((x == 0 || x == xSize-1 || y == 0 || y == ySize - 1 || z == 0 || z == zSize -1) && (scalarField[x][y][z] >= 0)){
                    out << std::fixed << setprecision(4) << -0.0005 << endl;
                }else{
				    out << std::fixed << setprecision(4) << scalarField[x][y][z] << " ";
                }
			}
            out << endl;
		}
        out << endl;
	}
    out.close();

	return 0;
}

void fitToGrid(float *point, bool isMin){
    //Get the sign
    int sign;
    if(isMin){
        sign = -1;
    }else{
        sign = 1;
    }
    for(int i = 0; i < 3; i++){
        if(point[i] < 0){       // Ensure the co-ordinate stays negative
            point[i] = -(abs(point[i]) + (sizeOfGrid - fmod(abs(point[i]), sizeOfGrid))) + (sign * 2 * sizeOfGrid);
        }  
        else{                   // Co-ordinate
            point[i] = point[i] + (sizeOfGrid - fmod(point[i], sizeOfGrid)) + (sign * 2 * sizeOfGrid);
        }
    }
}

float distanceToTriangle(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C, unsigned int faceID){
    // Check that the point does not lie on one of the vertices
    if(P == A || P == B || P == C){
        return 0;
    }

    // Can be pre-calculated
    glm::vec3 AB = B - A;
    glm::vec3 AC = C - A;
    glm::vec3 BA = A - B;
    glm::vec3 BC = C - B;
    glm::vec3 CA = A - C;
    glm::vec3 CB = B - C;

    // Calculate the normal of the triangle
    glm::vec3 N = fNormals[faceID];

    // Determine if inside the triangle
    glm::vec3 AP = P - A;
    glm::vec3 BP = P - B;
    glm::vec3 CP = P - C;
    float testAB = glm::dot(glm::cross(AB, N), AP);
    float testBC = glm::dot(glm::cross(BC, N), BP);
    float testCA = glm::dot(glm::cross(CA, N), CP);
    //cout << "Tests: " << testAB << " " << testBC << " " << testCA << endl;
    if(testAB <= 0 && testBC <= 0 && testCA <= 0){
        // Inside the triangle
        // Angle between normal and AP
        glm::vec3 AP = glm::normalize(P - A);
        float alpha = glm::dot(AP, N);

        // Vector P to P projected on the plane
        float dist2Plane = glm::length(A - P) * alpha;
        float s = sign(glm::dot(A - P, N));

        return s * abs(dist2Plane);
    }
    
    // Check for closest to vertex A
    float dAB = glm::dot(AB,AP);
    float dAC = glm::dot(AC,AP);
    if(dAB <= 0 && dAC <= 0){
        // Calculate the sign using vertex normal
        float s = sign(glm::dot(A - P, inputMesh.vNormals[inputMesh.faces[faceID].x]));
        float d = glm::length(A - P);

        return s * d;
    }

    // Check for closest to vertex B
    float dBA = glm::dot(BA,BP);
    float dBC = glm::dot(BC,BP);
    if(dBA <= 0 && dBC <= 0){
        // Calculate the sign using vertex normal
        float s = sign(glm::dot(B - P, inputMesh.vNormals[inputMesh.faces[faceID].y]));
        float d = glm::length(B - P);

        return s * d;
    }

    // Check for closest to vertex C
    float dCA = glm::dot(CA,CP);
    float dCB = glm::dot(CB,CP);
    if(dCA <= 0 && dCB <= 0){
        // Calculate the sign using vertex normal
        float s = sign(glm::dot(C - P, inputMesh.vNormals[inputMesh.faces[faceID].z]));
        float d = glm::length(C - P);

        return s * d;
    }
    
    // Check for closest to edge AB
    if(dAB >= 0 && dBA >= 0 && testAB >= 0){
        float d = distanceToEdge(P, A, B);
        float s = sign(glm::dot(A-P, eNormals[faceID * 3]));

        return s * d;
    }

    // Check for closest to edge BC
    if(dBC >= 0 && dCB >= 0 && testBC >= 0){
        float d = distanceToEdge(P, B, C);
        float s = sign(glm::dot(B-P, eNormals[faceID * 3 + 1]));
        return s * d;
    }

    // Check for closest to edge CA
    if(dCA >= 0 && dAC >= 0 && testCA >= 0){
        float d = distanceToEdge(P, C, A);
        float s = sign(glm::dot(C-P, eNormals[faceID * 3 + 2]));
        return s * d;
    }

    return 0;
}

float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B){
    return glm::length(glm::cross(P-A, P-B)) / glm::length(B-A);
}

int sign(float value){
    if(value < 0){
        return -1;
    }
    else if(value > 0){
        return 1;
    }
    else{
        return -1;
    }
}

glm::vec3 getVertexNormal(unsigned int vertexID){
    vector<glm::vec3> triangleNormals;
    // Find the triangles using the vertex
    for(size_t fID = 0; fID < inputMesh.faces.size(); fID++){
        if(inputMesh.faces[fID].x == vertexID || inputMesh.faces[fID].y == vertexID || inputMesh.faces[fID].z == vertexID){
            triangleNormals.push_back(fNormals[fID]);
        }
    }
    // Average the triangle normals
    glm::vec3 normal = glm::vec3(0,0,0);
    for(glm::vec3 n: triangleNormals){
        normal = normal + n;
    }
    normal = glm::normalize(normal / float(triangleNormals.size()));

    return normal;
}

glm::vec3 findOtherHalfNormal(unsigned int v0, unsigned int v1){
    // Look for v1 to v0
    for(size_t fID = 0; fID < inputMesh.faces.size(); fID++){
        if((inputMesh.faces[fID].x == v1 && inputMesh.faces[fID].y == v0) || 
            (inputMesh.faces[fID].y == v1 && inputMesh.faces[fID].z == v0) ||
            (inputMesh.faces[fID].z == v1 && inputMesh.faces[fID].x == v0)
        ){
            return fNormals[fID];
        }
    }
    return glm::vec3(0);
}