#include "Simplification.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define DEFAULT_STOPPING_ERROR 5.0

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

    // Get the optional inputs from the user
    int desiredTriangleCount = 0;
    float stoppingError = DEFAULT_STOPPING_ERROR;
    bool answered = false;
    string answer;
    while (!answered)
    {
        cout << "Would you like to set a desired triangle count? (y/n)";
        cin >> answer;
        if(answer == "y" || answer == "Y"){
            cout << "Enter an integer for the number of triangles.";
            cin >> desiredTriangleCount;
            answered = true;
        }
        if(answer == "n" || answer == "N"){
            answered = true;
        }
    }
    answered = false;
    while (!answered)
    {
        cout << "Would you like to set a stopping error? (Default is 5.0) (y/n)";
        cin >> answer;
        if(answer == "y" || answer == "Y"){
            cout << "Enter an floating point number for the stopping error.";
            cin >> stoppingError;
            answered = true;
        }
        if(answer == "n" || answer == "N"){
            answered = true;
        }
    }
    
    cout << "Number of triangles: " << mesh.edges.size() / 3 << endl;
    cout << "Number of vertices: " << mesh.vertices.size() / 3 << endl;

    // TIMER 1 START
    auto start1 = std::chrono::high_resolution_clock::now();

    /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////
    cout << "Calculating Q for each vertex" << endl;
    // Calculate the Quadric Error Metric for each half-edge
    // Find Q for each vertex
    quadrics.resize(mesh.vertices.size());
    for(size_t vID = 0; vID < mesh.vertices.size(); vID++){
        // Update Q
        updateQ(vID);
        // Find one ring
        oneRings.push_back(findOneRing(vID));
    }

    cout << "Calculating error for each edge" << endl;
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
    
    // Complete the collapse in the order of lowest error cost first
    while(!collapseOrder.empty()){
        // Stop if the triangle count has been met
        if(int(mesh.edges.size() / 3) <= desiredTriangleCount){
            cout << "Stopping due to triangle count being met" << endl;
            break;
        }
        // Stop if the error is too great
        if(collapseOrder[0].first > stoppingError && desiredTriangleCount == 0){
            cout << "Stopping due to error cost being too great" << endl;
            break;
        }
        // Get the edge to collapse
        unsigned int edge = collapseOrder[0].second;
        unsigned int otherEdge = mesh.otherhalves[edge];

        // Get the vertices involved
        unsigned int keptVertex = mesh.edges[edge]; 
        unsigned int goneVertex = mesh.edges[otherEdge];
        glm::vec3 newVertexPosition = optimalVertexPosition[collapseOrder[0].second];        

        // Remove the faces, swap the order if the second face is the second to last in the array
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
            if(oneRings[vID].find(goneVertex) != oneRings[vID].end()){
                oneRings[vID].erase(oneRings[vID].find(goneVertex));
            }
        }
        // Remove the old vertex

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
        for(auto id: oneRings[keptVertex]) {
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
        
        // Remove the error costs
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

        if(counter % 1000 == 0){
            cout << counter << " iterations completed." << endl;
        }
        counter++;
    }


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

    outputToDiredge();

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
    // TIMER 1 END
    auto end1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);

    cout << "Outputting. " << counter << " iterations completed. Time taken: " << float(duration1.count()) / 1000000.0 << " seconds." << endl;

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
    else{   // Assign -1 to the removed edges otherhalves
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

void updateQ(unsigned int vertexID){
    vector<glm::mat4> allK;
    glm::mat4 Q = glm::mat4(0);

    // Find the 1-ring
    for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
        if(mesh.edges[edgeID] == vertexID){
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
    // Check if it is a "valid" collapse
    int intersectionCount = 0;
    // Count one-ring intersections
    for(auto id: oneRings[mesh.edges[edgeID]]) {
        if(oneRings[mesh.edges[mesh.otherhalves[edgeID]]].find(id) != oneRings[mesh.edges[mesh.otherhalves[edgeID]]].end()){
            intersectionCount++;
        }
    }

    // If the one ring of both vertices intersect at more than 2 vertices it is an invalid operation
    if (intersectionCount > 2){
        // Invalid collapse operation
        return -1;
    }
    // If the object is a tetrahedron it is invalid to collapse 
    else if ((intersectionCount == 2) && (oneRings[mesh.edges[edgeID]].size() == 3) && (oneRings[mesh.edges[mesh.otherhalves[edgeID]]].size() == 3)) {
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
    // Check if invertible
    if(glm::determinant(newMat) == 0){
        // Place on the from vertex
        vertexPos = glm::vec4(mesh.vertices[mesh.edges[edgeID]], 1);
    }
    else{
        // Optimal vertex position
        vertexPos = glm::inverse(newMat)*glm::vec4(0,0,0,1);
    }

    // Check position does not equal NaN
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
    cout << "Outputting to lowpoly.diredge" << endl;
    // Write the output to a directed edge file format
	ofstream out("lowpoly.diredge");
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
    cout << "Outputting to lowpoly.obj" << endl;
    // Write to .obj file
    ofstream out("lowpoly.obj");
	for(size_t i = 0; i < mesh.vertices.size(); i++){
		out << "v " << std::fixed << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << endl;
	}
    for(size_t i = 0; i < mesh.vertexNormals.size(); i++){
		out << "vn " << std::fixed << mesh.vertexNormals[i][0] << " " << mesh.vertexNormals[i][1] << " " << mesh.vertexNormals[i][2] << endl;
	}
	for(size_t i = 0; i < mesh.edges.size(); i+=3){
		out << "f " 
        << mesh.edges[i] + 1 << "//" << mesh.edges[i] + 1 << " "
        << mesh.edges[i+1] + 1 << "//" << mesh.edges[i+1] + 1 << " "
        << mesh.edges[i+2] + 1 << "//" << mesh.edges[i+2] + 1 
        << endl;
	}
	out.close();
}

int checkForFlip(unsigned int goneVertexID, unsigned int keptVertexID, glm::vec3 newVertexPosition){
    volatile bool found = false;
    // Loop through all the faces to find each triangle the vertex is used in
    #pragma omp parallel for shared(found)
    for(size_t edgeID = 0; edgeID < mesh.edges.size(); edgeID++){
        if(found){
            continue;
        }

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
                    found = true;
                }
                glm::vec3 newNorm = glm::normalize(glm::cross(B-A, C-A));
                // Check if the normal has been flipped
                if(glm::dot(oldNorm, newNorm) < 0){
                    found = true;
                }
            }
        }
    }

    if(found){
        return -1;
    }
    return 0;
}