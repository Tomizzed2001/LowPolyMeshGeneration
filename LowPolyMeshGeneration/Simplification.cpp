#include "Simplification.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define DESIRED_TRIANGLE_COUNT 100

int main(int argc, char** argv){
    // Read in the .diredge file provided
    std::ifstream inFile(argv[1]);
    if(inFile.is_open()){
        string newLine;
        string type;
        glm::vec3 vec;
        int id;
        // Loop through all lines of the file
        while(getline(inFile, newLine)){
            // Convert line to string stream and get identifier
            stringstream line(newLine);
            line >> type;
            // Check if vertex line
            if(type == "v"){
                line >> vec.x >> vec.y >> vec.z;
                vertices.push_back(vec);
                continue;
            }
            // Check if first directed edge line
            if(type == "fde"){
                line >> id;
                firstDirectedEdges.push_back(id);
                continue;
            }
            // Check if face 
            if(type == "f"){
                line >> vec.x >> vec.y >> vec.z;
                faces.push_back(vec.x);
                faces.push_back(vec.y);
                faces.push_back(vec.z);
                continue;
            }
            // Check if other half
            if(type == "oh"){
                line >> id;
                otherHalf.push_back(id);
                continue;
            }
        }
    }
    inFile.close();

    // Calculate the Quadric Error Metric for each half-edge
    // Step 1: Find Q for each vertex
    quadrics.resize(vertices.size());
    for(size_t vID = 0; vID < vertices.size(); vID++){
        // Vector to store the K values of the 1-ring
        vector<glm::mat4> allK;

        // Find the 1-ring
        for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
            if(faces[edgeID] == vID){
                // Calculate the triangle ID and return K
                // allK.push_back(findK(vID / 3));
                break;
            }
        }
    }

    findK(4);



}

glm::mat4 findK(int triangleID){
    // Get the triangle
    glm::vec3 A, B, C;
    A = vertices[faces[triangleID * 3]];
    B = vertices[faces[triangleID * 3 + 1]];
    C = vertices[faces[triangleID * 3 + 2]];

    // Get the plane
    // Step 1: Calculate the normal
    glm::vec3 AB = B-A;
    glm::vec3 AC = C-A;
    glm::vec3 norm = glm::normalize(glm::cross(AB, AC));
    // Step 2: Calculate the offset
    float d = -(norm.x * A.x + norm.y * A.y + norm.z * A.z);
    // Step 3: Plane is [a,b,c,d]T
    glm::vec4 p = glm::vec4(norm.x, norm.y, norm.z, d);
    // Step 4: Build the matrix
    glm::mat4 K = glm::outerProduct(p,p);

    return glm::mat4(K);
}