#include "Simplification.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define DESIRED_TRIANGLE_COUNT 500

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
                mesh.vertices.push_back(vec);
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
                mesh.edges.push_back(vec.x);
                mesh.edges.push_back(vec.y);
                mesh.edges.push_back(vec.z);
                continue;
            }
            // Check if other half
            if(type == "oh"){
                line >> id;
                mesh.otherhalves.push_back(id);
                continue;
            }
        }
    }
    inFile.close();

    outputToObject(0);

    cout << "Number of triangles: " << mesh.edges.size() / 3 << endl;
    cout << "Number of mesh.vertices: " << mesh.vertices.size() / 3 << endl;

    /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////

    // Calculate the Quadric Error Metric for each half-edge
    // Find Q for each vertex
    quadrics.resize(mesh.vertices.size());
    for(size_t vID = 0; vID < mesh.vertices.size(); vID++){
        // Update Q
        updateQ(vID);
        // Find one ring
        oneRings.push_back(findOneRing(vID));
    }

    // Make a vector of the collapse order as a set of pairs <error, edgeID>
    optimalVertexPosition.resize(mesh.otherhalves.size());
    vector<pair<float, int>> collapseOrder;
    for(size_t ohID = 0; ohID < mesh.otherhalves.size(); ohID++){
        // Get the error cost for a valid half edge
        errorCosts.push_back(getEdgeError(ohID));
        if(errorCosts.back() != -1){
            collapseOrder.push_back(make_pair(errorCosts.back(), ohID));
        }
    }
    
    // Sort with smallest length first
    sort(collapseOrder.begin(), collapseOrder.end());

    /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////

    cout << "Beginning main loop" << endl;
    vector<unsigned int> removedVertices;
    int counter = 1;

    // TIMER 1 START
    auto start1 = std::chrono::high_resolution_clock::now();

    // Complete the collapse in the order of lowest error cost first
    while((mesh.edges.size() / 3) > DESIRED_TRIANGLE_COUNT && !collapseOrder.empty()){
        // Get the edge to collapse
        unsigned int edge = collapseOrder[0].second;
        unsigned int otherEdge = mesh.otherhalves[edge];

        // Get the vertices involved
        unsigned int keptVertex = mesh.edges[edge]; 
        unsigned int goneVertex = mesh.edges[otherEdge];
        glm::vec3 newVertexPosition = optimalVertexPosition[collapseOrder[0].second];        

        if ((otherEdge / 3) * 3 == mesh.edges.size()-3){
            removeFace(otherEdge / 3, 1);
            removeFace(edge / 3, 2);
        }
        else{
            removeFace(edge / 3, 1);
            removeFace(otherEdge / 3, 2);
        }

        // Update the one rings
        for(auto vID: oneRings[goneVertex]) {
            if((oneRings[keptVertex].find(vID) == oneRings[keptVertex].end()) && (vID != keptVertex)){
                oneRings[keptVertex].insert(vID);
                oneRings[vID].insert(keptVertex);
            }
        }
        // Remove the old vertex
        oneRings[keptVertex].erase(oneRings[keptVertex].find(goneVertex));

        vector<int> firstEdge, secondEdge;

        // Replace all instances of the goneVertex with the keptVertex
        for(size_t eID = 0; eID < mesh.edges.size(); eID++){
            // Check for removed vertex
            if(mesh.edges[eID] == goneVertex){
                mesh.edges[eID] = keptVertex;
                updatedEdges.insert(eID);
                updatedEdges.insert(mesh.otherhalves[eID]);
            }
            else if(mesh.edges[eID] == keptVertex){
                updatedEdges.insert(eID);
                updatedEdges.insert(mesh.otherhalves[eID]);
            }
            // Check for removed half edge on the first removed face
            if(mesh.otherhalves[eID] == -1){
                firstEdge.push_back(eID);
            }
            else if(mesh.otherhalves[eID] == -2){
                secondEdge.push_back(eID);
            }
        }

        // Update the other halves for the kept edges of the removed triangles
        mesh.otherhalves[firstEdge[0]] = firstEdge[1];
        mesh.otherhalves[firstEdge[1]] = firstEdge[0];
        mesh.otherhalves[secondEdge[0]] = secondEdge[1];
        mesh.otherhalves[secondEdge[1]] = secondEdge[0];

        // Add removed vertex to an array so the hole set can be updated at the end
        removedVertices.push_back(goneVertex);

        // Update the position of the vertex to the optimal one
        mesh.vertices[keptVertex] = newVertexPosition;

        // Update the Quadric for the kept vertex
        quadrics[keptVertex] = quadrics[keptVertex] + quadrics[goneVertex];

        // Get all the edges in the 2-ring of the vertex to update them
        std::unordered_set<unsigned int> oneRing = findOneRing(keptVertex);
        for(auto id: oneRing) {
            for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
                if(mesh.edges[edgeID] == id){
                    // Vertex is A in the triangle
                    if(edgeID % 3 == 0){
                        updatedEdges.insert(edgeID);
                        updatedEdges.insert(edgeID + 2);
                        continue;
                    }
                    else{
                        updatedEdges.insert(edgeID);
                        updatedEdges.insert(edgeID - 1);
                    }
                }
            }
        }

        
        // Update the cost for all edges involving the kept vertex
        for(auto edge: updatedEdges){
            if(edge >= 0){
                errorCosts[edge] = getEdgeError(edge);
            }
        }


        updatedEdges.clear();

        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();

        // Re-Generate Collapse Order
        collapseOrder.clear();
        for(size_t ohID = 0; ohID < mesh.otherhalves.size(); ohID++){
            if(errorCosts[ohID] != -1){
                collapseOrder.push_back(make_pair(errorCosts[ohID], ohID));
            }
        }

        // Sort with smallest length first
        sort(collapseOrder.begin(), collapseOrder.end());

        // DEBUG OBJECT FILE
        if(counter % 100 == 0){
            // TIMER 1 END
            auto end1 = std::chrono::high_resolution_clock::now();

            auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

            cout << "Outputting. " << counter << " iterations completed. Time taken: " << float(duration1.count()) / 1000000.0 << " seconds." << endl;
            outputToObject(counter);

            start1 = std::chrono::high_resolution_clock::now();
        }
        counter++;
    }

    outputToDiredge();

    cout << "Removing vertices" << endl;
    // Remove un-used vertices and cascade the changes
    for(size_t i = 0; i < removedVertices.size(); i++){
        // Remove the vertex from vertices
        mesh.vertices.erase(mesh.vertices.begin() + (removedVertices[i]));
        // Decrement all vertices that had higher index
        for(size_t eID = 0; eID < mesh.edges.size(); eID++){
            if(mesh.edges[eID] > removedVertices[i]){
                mesh.edges[eID]--;
            }
        }
        // Decrement the removed vertex list
        for(size_t j = i; j < removedVertices.size(); j++){
            if(removedVertices[j] > removedVertices[i]){
                removedVertices[j]--;
            }
        }
    }

    cout << "Producing Normals" << endl;
    // Produce the smooth vertex normals for the output
    for(size_t vID = 0; vID < mesh.vertices.size(); vID++){
        vector<glm::vec3> triangleNormals;
        // Find the triangles using the vertex
        for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
            if(mesh.edges[edgeID] == vID){
                int triangleID = edgeID / 3;
                // Get the triangle
                glm::vec3 A, B, C;
                A = mesh.vertices[mesh.edges[triangleID * 3]];
                B = mesh.vertices[mesh.edges[triangleID * 3 + 1]];
                C = mesh.vertices[mesh.edges[triangleID * 3 + 2]];

                // Calculate the normal
                glm::vec3 norm = glm::normalize( glm::cross(B-A, C-A) );
                //cout << "New Normal: " << glm::to_string(norm) << endl;
                triangleNormals.push_back(norm);
            }
        }
        // Average the triangle normals
        glm::vec3 normal = glm::vec3(0,0,0);
        for(glm::vec3 n: triangleNormals){
            normal = normal + n;
        }
        normal = glm::normalize(normal / float(triangleNormals.size()));
        mesh.vertexNormals.push_back(normal);
    }

    cout << "Complete" << endl;

    // Output to .obj file
    outputToObject();
}

glm::mat4 findK(int triangleID){
    // Get the triangle
    glm::vec3 A = mesh.vertices[mesh.edges[triangleID * 3]];
    glm::vec3 B = mesh.vertices[mesh.edges[triangleID * 3 + 1]];
    glm::vec3 C = mesh.vertices[mesh.edges[triangleID * 3 + 2]];
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

std::unordered_set<unsigned int> findOneRing(unsigned int vertexID){
    std::unordered_set<unsigned int> oneRing;

    // Loop through the mesh.edges finding each triangle the vertex is part of
    for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
        if(mesh.edges[edgeID] == vertexID){
            // Vertex is A in the triangle
            if(edgeID % 3 == 0){
                oneRing.insert(mesh.edges[edgeID + 1]);
                oneRing.insert(mesh.edges[edgeID + 2]);
                continue;
            }
            // Vertex is B in the triangle
            if(edgeID % 3 == 1){
                oneRing.insert(mesh.edges[edgeID - 1]);
                oneRing.insert(mesh.edges[edgeID + 1]);
                continue;
            }
            // Vertex is C in the triangle
            if(edgeID % 3 == 2){
                oneRing.insert(mesh.edges[edgeID - 2]);
                oneRing.insert(mesh.edges[edgeID - 1]);
                continue;
            }
        }
    }

    return oneRing;
}

void removeFace(unsigned int triangleID, int faceNum){
    int faceID = triangleID * 3;
    //cout << "FaceID: " << faceID << endl;
    // Check that its not the last triangle in the array
    if (faceID != int(mesh.edges.size()-3)){
        // Swap the faces
        mesh.edges[faceID] = mesh.edges[mesh.edges.size()-3];
        mesh.edges[faceID+1] = mesh.edges[mesh.edges.size()-2];
        mesh.edges[faceID+2] = mesh.edges[mesh.edges.size()-1];
        updatedEdges.insert(faceID);
        updatedEdges.insert(faceID+1);
        updatedEdges.insert(faceID+2);

        // For each edge in the triangle
        for(int i = 0; i < 3; i++){
            // If the edge will be the same but flipped remove the edge
            if(faceID+i == mesh.otherhalves[mesh.otherhalves.size()-(3-i)]){
                mesh.otherhalves[faceID+i] = -1 * faceNum;
            }
            else{
                // If the edge already has no other half
                if(mesh.otherhalves[faceID+i] != -1){
                    mesh.otherhalves[mesh.otherhalves[faceID+i]] = -1 * faceNum;
                }
                // Swap the other half
                mesh.otherhalves[faceID+i] = mesh.otherhalves[mesh.otherhalves.size()-(3-i)];
                // Update the corresponding other half with the new position
                mesh.otherhalves[mesh.otherhalves[faceID+i]] = faceID+i;
            }
        }

    }
    else{
        mesh.otherhalves[mesh.otherhalves[mesh.edges.size()-3]] = -1 * faceNum;
        mesh.otherhalves[mesh.otherhalves[mesh.edges.size()-2]] = -1 * faceNum;
        mesh.otherhalves[mesh.otherhalves[mesh.edges.size()-1]] = -1 * faceNum;
    }

    // Remove the final face
    mesh.edges.pop_back();
    mesh.edges.pop_back();
    mesh.edges.pop_back();
    mesh.otherhalves.pop_back();
    mesh.otherhalves.pop_back();
    mesh.otherhalves.pop_back();
}

void findOtherHalf(unsigned int edgeID){
    unsigned int v0, v1;
    if (edgeID % 3 == 2){
        v0 = mesh.edges[edgeID];
        v1 = mesh.edges[edgeID-2];
    }
    else{
        v0 = mesh.edges[edgeID];
        v1 = mesh.edges[edgeID+1];
    }

    // Look for v1 to v0
    for(size_t fID = 0; fID < mesh.edges.size(); fID+=3){
        if(mesh.edges[fID] == v1 && mesh.edges[fID+1] == v0){
            mesh.otherhalves[edgeID] = fID + 0;
            break;
        }
        else if(mesh.edges[fID+1] == v1 && mesh.edges[fID+2] == v0){
            mesh.otherhalves[edgeID] = fID + 1;
            break;
        }
        else if(mesh.edges[fID+2] == v1 && mesh.edges[fID] == v0){
            mesh.otherhalves[edgeID] = fID + 2;
            break;
        }
    }

    mesh.otherhalves[mesh.otherhalves[edgeID]] = edgeID;

    return;
}

void updateQ(unsigned int vertexID){
    //cout << "Getting Q for: " << vertexID << endl;
    vector<glm::mat4> allK;
    glm::mat4 Q = glm::mat4(0);

    // Find the 1-ring
    for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
        if(mesh.edges[edgeID] == vertexID){
            //cout << "Triangle in one-ring: " << edgeID / 3 << endl;
            // Calculate the triangle ID and return K
            allK.push_back(findK(edgeID / 3));
        }
    }

    // Sum the K matrices to get Q for the vertex
    for(size_t i = 0; i < allK.size(); i++){
        Q = Q + allK[i];
    }

    quadrics[vertexID] = Q;

    return;
}

float getEdgeError(unsigned int edgeID){
    // cout << edgeID << endl;
    // Check if it is a "valid" collapse
    std::unordered_set<unsigned int> aOneRing = findOneRing(mesh.edges[edgeID]);
    std::unordered_set<unsigned int> bOneRing = findOneRing(mesh.edges[mesh.otherhalves[edgeID]]);
    int intersectionCount = 0;
    for(auto id: aOneRing) {
        if(bOneRing.find(id) != bOneRing.end()){
            intersectionCount++;
        }
    }
    // If the one ring of both vertices intersect at more than 2 vertices it is an invalid operation
    if (intersectionCount > 2){
        // Invalid collapse operation
        return -1;
    }

    // Get the quadric for each vertex
    glm::mat4 Q1 = quadrics[mesh.edges[edgeID]];
    glm::mat4 Q2 = quadrics[mesh.edges[mesh.otherhalves[edgeID]]];

    // Get optimal vertex position for this collapse
    glm::mat4 newMat = Q1+Q2;
    newMat[0][3] = 0;
    newMat[1][3] = 0;
    newMat[2][3] = 0;
    newMat[3][3] = 1;

    glm::vec4 vertexPos;
    if(glm::determinant(newMat) == 0){
        //cout << "Not Invertible" << endl;
        // Place on the from vertex
        vertexPos = glm::vec4(mesh.vertices[mesh.edges[edgeID]], 1);
    }
    else{
        vertexPos = glm::inverse(newMat)*glm::vec4(0,0,0,1);
    }

    if(vertexPos.x != vertexPos.x){
        vertexPos = glm::vec4(mesh.vertices[mesh.edges[edgeID]], 1);
    }

    optimalVertexPosition[edgeID] = glm::vec3(vertexPos);


    // Using the new vertex determine if this will cause any triangle flips (self-intersections)
    if(checkForFlip(mesh.edges[mesh.otherhalves[edgeID]] , mesh.edges[edgeID], optimalVertexPosition[edgeID]) == -1
        || checkForFlip( mesh.edges[edgeID], mesh.edges[mesh.otherhalves[edgeID]] ,optimalVertexPosition[edgeID]) == -1){
        return -1;
    }

    return abs(glm::dot(vertexPos, ((Q1 + Q2) * vertexPos)));
}

void outputToDiredge(){
    cout << "Dumping to dump.diredge" << endl;
    // Write the output to a directed edge file format
	ofstream out("dump.diredge");
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
	for(size_t i = 0; i < firstDirectedEdges.size(); i++){
		out << "fde " << firstDirectedEdges[i] << endl;
	} 
	for(size_t i = 0; i < mesh.edges.size(); i+=3){
		out << "f " << mesh.edges[i] << " " << mesh.edges[i+1] << " " << mesh.edges[i+2] << endl;
	}
	for(size_t i = 0; i < mesh.otherhalves.size(); i++){
		out << "oh " << mesh.otherhalves[i] << endl;
	}

	out.close();
}

void outputToObject(){
    ofstream out("out.obj");
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
    for(size_t i = 0; i < mesh.vertexNormals.size(); i++){
		out << "vn " << std::fixed << mesh.vertexNormals[i][0] << " " << mesh.vertexNormals[i][1] << " " << mesh.vertexNormals[i][2] << endl;
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

void outputToObject(int num){
    std::string var = "outputs/out" + to_string(num) + ".obj";
    ofstream out(var);
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
    for(size_t i = 0; i < mesh.vertexNormals.size(); i++){
		out << "vn " << std::fixed << mesh.vertexNormals[i][0] << " " << mesh.vertexNormals[i][1] << " " << mesh.vertexNormals[i][2] << endl;
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

float getEdgeLength(unsigned int edgeID){
    // Check if it is a "valid" collapse
    std::unordered_set<unsigned int> aOneRing = findOneRing(mesh.edges[edgeID]);
    std::unordered_set<unsigned int> bOneRing = findOneRing(mesh.edges[mesh.otherhalves[edgeID]]);
    int intersectionCount = 0;
    for(auto id: aOneRing) {
        if(bOneRing.find(id) != bOneRing.end()){
            intersectionCount++;
        }
    }
    // If the one ring of both vertices intersect at more than 2 vertices it is an invalid operation
    if (intersectionCount > 2){
        // Invalid collapse operation
        return -1;
    }

    // Get the vertices
    glm::vec3 A = mesh.vertices[mesh.edges[edgeID]];
    glm::vec3 B;
    if(edgeID % 3 == 2){
        B = mesh.vertices[mesh.edges[edgeID-2]];
    }
    else{
        B = mesh.vertices[mesh.edges[edgeID+1]];
    }

    return glm::length(B-A);
}

int checkForFlip(unsigned int goneVertexID, unsigned int keptVertexID, glm::vec3 newVertexPosition){
    // Loop through all the faces to find each triangle the vertex is used in
    for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
        // Find a triangle
        if(mesh.edges[edgeID] == goneVertexID){
            int triangleID = (edgeID / 3) * 3;
            glm::vec3 A = mesh.vertices[mesh.edges[triangleID]];
            glm::vec3 B = mesh.vertices[mesh.edges[triangleID + 1]];
            glm::vec3 C = mesh.vertices[mesh.edges[triangleID + 2]];
            // If the triangle will be removed skip this one
            if(mesh.edges[triangleID] == keptVertexID || mesh.edges[triangleID + 1] == keptVertexID || mesh.edges[triangleID + 2] == keptVertexID){
                continue;
            }
            else{
                // Get the normal of the triangle
                glm::vec3 oldNorm = glm::normalize(glm::cross(B-A, C-A));
                if(mesh.edges[triangleID] == goneVertexID){
                    A = newVertexPosition;
                }
                else if(mesh.edges[triangleID + 1] == goneVertexID){
                    B = newVertexPosition;
                }
                else{
                    C = newVertexPosition;
                }
                if(glm::cross(B-A, C-A) == glm::vec3(0,0,0)){
                    return -1;
                }
                glm::vec3 newNorm = glm::normalize(glm::cross(B-A, C-A));
                // Check if the normal has been flipped
                if(glm::dot(oldNorm, newNorm) < 0){
                    return -1;
                }
            }
        }
    }

    return 0;
}