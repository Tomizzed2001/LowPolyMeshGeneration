#include "Simplification.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define DESIRED_TRIANGLE_COUNT 6000

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

    outputToObject(0);

    cout << "Number of triangles: " << faces.size() / 3 << endl;
    cout << "Number of vertices: " << vertices.size() / 3 << endl;

    ///////////////////////////////////////////// TEST THE EDGE COLLAPSE ////////////////////////////////////////////
    /*
    // Get the edge to collapse
    unsigned int edge = 0;
    unsigned int otherEdge = otherHalf[edge];

    //cout << "Edge: " << edge << " Other Edge: " << otherEdge << endl;

    // Get the vertices involved
    unsigned int keptVertex = faces[edge]; 
    unsigned int goneVertex = faces[otherEdge];

    //cout << "Keep Vertex: " << keptVertex << " Remove Vertex: " << goneVertex << endl;

    if ((otherEdge / 3) * 3 == faces.size()-3){
        cout << "OTHER IS ON THE FINAL FACE. First Face: " << otherEdge / 3 << " Second Face: " << edge / 3 << endl;
        removeFaceNew(otherEdge / 3, 1);
        removeFaceNew(edge / 3, 2);
    }
    else{
        cout << "First Face: " << edge / 3 << " Second Face: " << otherEdge / 3 << endl;
        removeFaceNew(edge / 3, 1);
        removeFaceNew(otherEdge / 3, 2);
    }

    vector<int> firstEdge, secondEdge; 

    // Replace all instances of the goneVertex with the keptVertex
    for(size_t eID = 0; eID < faces.size(); eID++){
        // Check for removed vertex
        if(faces[eID] == goneVertex){
            faces[eID] = keptVertex;
        }
        // Check for removed half edge on the first removed face
        if(otherHalf[eID] == -1){
            firstEdge.push_back(eID);
        }
        else if(otherHalf[eID] == -2){
            secondEdge.push_back(eID);
        }
    }

    otherHalf[firstEdge[0]] = firstEdge[1];
    otherHalf[firstEdge[1]] = firstEdge[0];
    otherHalf[secondEdge[0]] = secondEdge[1];
    otherHalf[secondEdge[1]] = secondEdge[0];

    removedVertices.push_back(goneVertex);

    outputToObject(1);

    outputToDiredge();

    return 0;
    */
    ///////////////////////////////////////////// TEST THE EDGE COLLAPSE ////////////////////////////////////////////

    /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////

    // Calculate the Quadric Error Metric for each half-edge
    // Find Q for each vertex
    quadrics.resize(vertices.size());
    for(size_t vID = 0; vID < vertices.size(); vID++){
        updateQ(vID);
    }

    // Make a vector of the collapse order as a set of pairs <error, edgeID>
    optimalVertexPosition.resize(otherHalf.size());
    vector<pair<float, int>> collapseOrder;
    for(size_t ohID = 0; ohID < otherHalf.size(); ohID++){
        // Get the error cost for a valid half edge
        errorCosts.push_back(getEdgeError(ohID));
        if(errorCosts.back() != -1){
            collapseOrder.push_back(make_pair(errorCosts.back(), ohID));
        }
    }

    /*
    for(size_t i = 0; i < errorCosts.size(); i++){
        cout << "Error of " << i << ": " << errorCosts[i] << endl;
    }
    */
    
    // Sort with smallest length first
    sort(collapseOrder.begin(), collapseOrder.end());
    
    /*
    for(size_t i = 0; i < collapseOrder.size(); i++){
        cout << collapseOrder[i].first << " " << collapseOrder[i].second << endl;
    }
    */

    /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////

    cout << "Beginning main loop" << endl;
    vector<unsigned int> removedVertices;
    int counter = 1;
    //cout << "Collapsing edges" << endl;
    // Complete the collapse in the order of lowest error cost first
    while((faces.size() / 3) > DESIRED_TRIANGLE_COUNT && !collapseOrder.empty()){
        //cout << "Run: " << counter << endl;
        // Get the edge to collapse
        unsigned int edge = collapseOrder[0].second;
        unsigned int otherEdge = otherHalf[edge];

        //cout << "Run: " << counter << " Edge: " << collapseOrder[0].second;

        // Get the vertices involved
        unsigned int keptVertex = faces[edge]; 
        unsigned int goneVertex = faces[otherEdge];
        glm::vec3 newVertexPosition = optimalVertexPosition[collapseOrder[0].second];
        //cout << " Kept Vertex: " << keptVertex << " Gone Vertex: " << goneVertex << " Optimal vertex: " << glm::to_string(newVertexPosition) << endl;

        

        if ((otherEdge / 3) * 3 == faces.size()-3){
            removeFaceNew(otherEdge / 3, 1);
            removeFaceNew(edge / 3, 2);
        }
        else{
            removeFaceNew(edge / 3, 1);
            removeFaceNew(otherEdge / 3, 2);
        }
 

        vector<int> firstEdge, secondEdge;

        // Replace all instances of the goneVertex with the keptVertex
        for(size_t eID = 0; eID < faces.size(); eID++){
            // Check for removed vertex
            if(faces[eID] == goneVertex){
                faces[eID] = keptVertex;
                updatedEdges.push_back(eID);
            }
            else if(faces[eID] == keptVertex){
                updatedEdges.push_back(eID);
            }
            // Check for removed half edge on the first removed face
            if(otherHalf[eID] == -1){
                firstEdge.push_back(eID);
            }
            else if(otherHalf[eID] == -2){
                secondEdge.push_back(eID);
            }
        }

        // Update the other halves for the kept edges of the removed triangles
        otherHalf[firstEdge[0]] = firstEdge[1];
        otherHalf[firstEdge[1]] = firstEdge[0];
        otherHalf[secondEdge[0]] = secondEdge[1];
        otherHalf[secondEdge[1]] = secondEdge[0];

        // Add removed vertex to an array so the hole set can be updated at the end
        removedVertices.push_back(goneVertex);

        // Update the position of the vertex to the optimal one
        //cout << "New position: " << glm::to_string(newVertexPosition) << endl;
        vertices[keptVertex] = newVertexPosition;

        /////////////////////////////////////////////// COLLAPSE ORDERING ///////////////////////////////////////////////

        // Update the Quadric for the kept vertex
        quadrics[keptVertex] = quadrics[keptVertex] + quadrics[goneVertex];

        // Get all the edges in the 2-ring of the vertex to update them
        std::unordered_set<unsigned int> oneRing = findOneRing(keptVertex);
        for(auto id: oneRing) {
            for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
                if(faces[edgeID] == id){
                    // Vertex is A in the triangle
                    if(edgeID % 3 == 0){
                        updatedEdges.push_back(edgeID);
                        updatedEdges.push_back(edgeID + 2);
                        continue;
                    }
                    else{
                        updatedEdges.push_back(edgeID);
                        updatedEdges.push_back(edgeID - 1);
                    }
                }
            }
        }

        // TIMER 1 START
        auto start1 = std::chrono::high_resolution_clock::now();
        
        // Update the cost for all edges involving the kept vertex
        for(size_t i = 0; i < updatedEdges.size(); i++){
            errorCosts[updatedEdges[i]] = getEdgeError(updatedEdges[i]);
            if(i >= 6){
                errorCosts[otherHalf[updatedEdges[i]]] = getEdgeError(otherHalf[updatedEdges[i]]);
            }
        }

        updatedEdges.clear();
        
        //cout << "Re-building" << endl;

        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();
        errorCosts.pop_back();

        /*
        cout << "Brute force list " << endl;
        for(size_t i = 0; i < errorCosts.size(); i++){
            float newCost = getEdgeError(i);
            if(errorCosts[i] != newCost){
                cout << i << " Vertex A: " << faces[i] << " Vertex B: " << faces[otherHalf[i]] << " Old cost: " << errorCosts[i] << " New cost: " << newCost << endl;
            }
            errorCosts[i] = newCost;
        }
        */

        // TIMER 1 END
        auto end1 = std::chrono::high_resolution_clock::now();

        auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
 
        //cout << "Time taken by function 1: " << duration1.count() << " microseconds" << endl;


        // Re-Generate Collapse Order
        collapseOrder.clear();
        for(size_t ohID = 0; ohID < otherHalf.size(); ohID++){
            if(errorCosts[ohID] != -1){
                collapseOrder.push_back(make_pair(errorCosts[ohID], ohID));
            }
        }

        // Sort with smallest length first
        sort(collapseOrder.begin(), collapseOrder.end());

        // DEBUG OBJECT FILE
        if(counter % 50 == 0){
            cout << "Outputting. " << counter << " iterations completed" << endl;
            outputToObject(counter);
        }
        counter++;
    }

    outputToDiredge();

    //outputToObject();
    //cout << "OUTPUT" << endl;

    cout << "Removing vertices" << endl;
    // Remove un-used vertices and cascade the changes
    for(size_t i = 0; i < removedVertices.size(); i++){
        // Remove the vertex from vertices
        vertices.erase(vertices.begin() + (removedVertices[i]));
        // Decrement all vertices that had higher index
        for(size_t eID = 0; eID < faces.size(); eID++){
            if(faces[eID] > removedVertices[i]){
                faces[eID]--;
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
    for(size_t vID = 0; vID < vertices.size(); vID++){
        vector<glm::vec3> triangleNormals;
        // Find the triangles using the vertex
        for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
            if(faces[edgeID] == vID){
                int triangleID = edgeID / 3;
                // Get the triangle
                glm::vec3 A, B, C;
                A = vertices[faces[triangleID * 3]];
                B = vertices[faces[triangleID * 3 + 1]];
                C = vertices[faces[triangleID * 3 + 2]];

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
        vertexNormals.push_back(normal);
        //cout << "Final Normal: " << glm::to_string(normal) << endl;
        //break;
    }

    cout << "Complete" << endl;
    // Output to .obj file
    outputToObject();
}

glm::mat4 findK(int triangleID){
    // Get the triangle
    glm::vec3 A = vertices[faces[triangleID * 3]];
    glm::vec3 B = vertices[faces[triangleID * 3 + 1]];
    glm::vec3 C = vertices[faces[triangleID * 3 + 2]];
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

    // Loop through the faces finding each triangle the vertex is part of
    for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
        if(faces[edgeID] == vertexID){
            // Vertex is A in the triangle
            if(edgeID % 3 == 0){
                oneRing.insert(faces[edgeID + 1]);
                oneRing.insert(faces[edgeID + 2]);
                continue;
            }
            // Vertex is B in the triangle
            if(edgeID % 3 == 1){
                oneRing.insert(faces[edgeID - 1]);
                oneRing.insert(faces[edgeID + 1]);
                continue;
            }
            // Vertex is C in the triangle
            if(edgeID % 3 == 2){
                oneRing.insert(faces[edgeID - 2]);
                oneRing.insert(faces[edgeID - 1]);
                continue;
            }
        }
    }

    return oneRing;
}

void removeFace(unsigned int triangleID){
    unsigned int faceID = triangleID * 3;
    // Check that its not the last triangle in the array and swap if not
    if (faceID != faces.size()-3){
        cout << "NOT LAST TRIANGLE" << endl;
        vector<pair<unsigned int, unsigned int>> swapMap;
        // Swap vertices
        swapMap.push_back(make_pair(faces[faceID], faces[faces.size()-3]));
        faces[faceID] = faces[faces.size()-3];
        swapMap.push_back(make_pair(faces[faceID+1], faces[faces.size()-2]));
        faces[faceID+1] = faces[faces.size()-2];
        swapMap.push_back(make_pair(faces[faceID+2], faces[faces.size()-1]));
        faces[faceID+2] = faces[faces.size()-1];

        cout << "OH: " << otherHalf[faces.size()-3] << endl;
        cout << "OH: " << otherHalf[faces.size()-2] << endl;
        cout << "OH: " << otherHalf[faces.size()-1] << endl;

        outputToDiredge();

        // Update otherHalf to match swapped faces
        for(pair<unsigned int, unsigned int> swap: swapMap){
            cout << swap.first << " " << swap.second << endl;
            int pairsFound = 0;
            unsigned int firstIndex, secondIndex;
            for(size_t i = 0; i < otherHalf.size() || pairsFound == 2; i++){
                if(swap.first == otherHalf[i]){
                    firstIndex = i;
                    // Check the edge isnt the same but swapped
                    if(firstIndex == swap.second){
                        break;
                    }
                    pairsFound++;
                    continue;
                }
                else if(swap.second == otherHalf[i]){
                    secondIndex = i;
                    // Check the edge isnt the same but swapped
                    if(secondIndex == swap.first){
                        //cout << "CHECK" << endl;
                        break;
                    }
                    pairsFound++;
                    continue;
                }
            }
            // Must already exist and was forced to break out of loop
            if(pairsFound != 2){
                continue;
            }
            // Update values
            otherHalf[swap.second] = firstIndex;
            otherHalf[swap.first] = secondIndex;
            otherHalf[firstIndex] = swap.second;
            otherHalf[secondIndex] = swap.first;
        }
    }

    // Remove the final face
    faces.pop_back();
    faces.pop_back();
    faces.pop_back();

    return;
}

void removeFaceNew(unsigned int triangleID, int faceNum){
    int faceID = triangleID * 3;
    //cout << "FaceID: " << faceID << endl;
    // Check that its not the last triangle in the array
    if (faceID != int(faces.size()-3)){
        // Swap the faces
        faces[faceID] = faces[faces.size()-3];
        faces[faceID+1] = faces[faces.size()-2];
        faces[faceID+2] = faces[faces.size()-1];
        updatedEdges.push_back(faceID);
        updatedEdges.push_back(faceID+1);
        updatedEdges.push_back(faceID+2);

        // For each edge in the triangle
        for(int i = 0; i < 3; i++){
            // If the edge will be the same but flipped remove the edge
            if(faceID+i == otherHalf[otherHalf.size()-(3-i)]){
                otherHalf[faceID+i] = -1 * faceNum;
            }
            else{
                // If the edge already has no other half
                if(otherHalf[faceID+i] != -1){
                    otherHalf[otherHalf[faceID+i]] = -1 * faceNum;
                }
                // Swap the other half
                otherHalf[faceID+i] = otherHalf[otherHalf.size()-(3-i)];
                // Update the corresponding other half with the new position
                otherHalf[otherHalf[faceID+i]] = faceID+i;
            }
        }

    }
    else{
        otherHalf[otherHalf[faces.size()-3]] = -1 * faceNum;
        otherHalf[otherHalf[faces.size()-2]] = -1 * faceNum;
        otherHalf[otherHalf[faces.size()-1]] = -1 * faceNum;
    }

    // Remove the final face
    faces.pop_back();
    faces.pop_back();
    faces.pop_back();
    otherHalf.pop_back();
    otherHalf.pop_back();
    otherHalf.pop_back();
}

void findOtherHalf(unsigned int edgeID){
    unsigned int v0, v1;
    if (edgeID % 3 == 2){
        v0 = faces[edgeID];
        v1 = faces[edgeID-2];
    }
    else{
        v0 = faces[edgeID];
        v1 = faces[edgeID+1];
    }

    // Look for v1 to v0
    for(size_t fID = 0; fID < faces.size(); fID+=3){
        if(faces[fID] == v1 && faces[fID+1] == v0){
            otherHalf[edgeID] = fID + 0;
            break;
        }
        else if(faces[fID+1] == v1 && faces[fID+2] == v0){
            otherHalf[edgeID] = fID + 1;
            break;
        }
        else if(faces[fID+2] == v1 && faces[fID] == v0){
            otherHalf[edgeID] = fID + 2;
            break;
        }
    }

    otherHalf[otherHalf[edgeID]] = edgeID;

    return;
}

void updateQ(unsigned int vertexID){
    //cout << "Getting Q for: " << vertexID << endl;
    vector<glm::mat4> allK;
    glm::mat4 Q = glm::mat4(0);

    // Find the 1-ring
    for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
        if(faces[edgeID] == vertexID){
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
    std::unordered_set<unsigned int> aOneRing = findOneRing(faces[edgeID]);
    std::unordered_set<unsigned int> bOneRing = findOneRing(faces[otherHalf[edgeID]]);
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
    glm::mat4 Q1 = quadrics[faces[edgeID]];
    glm::mat4 Q2 = quadrics[faces[otherHalf[edgeID]]];

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
        vertexPos = glm::vec4(vertices[faces[edgeID]], 1);
    }
    else{
        vertexPos = glm::inverse(newMat)*glm::vec4(0,0,0,1);
    }

    if(vertexPos.x != vertexPos.x){
        vertexPos = glm::vec4(vertices[faces[edgeID]], 1);
    }

    optimalVertexPosition[edgeID] = glm::vec3(vertexPos);


    // Using the new vertex determine if this will cause any triangle flips (self-intersections)
    if(checkForFlip(faces[otherHalf[edgeID]] , faces[edgeID], optimalVertexPosition[edgeID]) == -1
        || checkForFlip( faces[edgeID], faces[otherHalf[edgeID]] ,optimalVertexPosition[edgeID]) == -1){
        return -1;
    }

    return abs(glm::dot(vertexPos, ((Q1 + Q2) * vertexPos)));
}

void outputToDiredge(){
    cout << "Dumping to dump.diredge" << endl;
    // Write the output to a directed edge file format
	ofstream out("dump.diredge");
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
}

void outputToObject(){
    ofstream out("out.obj");
	for(size_t i = 0; i < vertices.size(); i++){
		out << "v " << std::fixed << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl;
	}
    for(size_t i = 0; i < vertexNormals.size(); i++){
		out << "vn " << std::fixed << vertexNormals[i][0] << " " << vertexNormals[i][1] << " " << vertexNormals[i][2] << endl;
	}
	for(size_t i = 0; i < faces.size(); i+=3){
		out << "f " 
        << faces[i] + 1 << "//" << " "
        << faces[i+1] + 1 << "//" << " "
        << faces[i+2] + 1 << "//"
        << endl;
	}
	out.close();
}

void outputToObject(int num){
    std::string var = "outputs/out" + to_string(num) + ".obj";
    ofstream out(var);
	for(size_t i = 0; i < vertices.size(); i++){
		out << "v " << std::fixed << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << endl;
	}
    for(size_t i = 0; i < vertexNormals.size(); i++){
		out << "vn " << std::fixed << vertexNormals[i][0] << " " << vertexNormals[i][1] << " " << vertexNormals[i][2] << endl;
	}
	for(size_t i = 0; i < faces.size(); i+=3){
		out << "f " 
        << faces[i] + 1 << "//" << " "
        << faces[i+1] + 1 << "//" << " "
        << faces[i+2] + 1 << "//"
        << endl;
	}
	out.close();
}

float getEdgeLength(unsigned int edgeID){
    // Check if it is a "valid" collapse
    std::unordered_set<unsigned int> aOneRing = findOneRing(faces[edgeID]);
    std::unordered_set<unsigned int> bOneRing = findOneRing(faces[otherHalf[edgeID]]);
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
    glm::vec3 A = vertices[faces[edgeID]];
    glm::vec3 B;
    if(edgeID % 3 == 2){
        B = vertices[faces[edgeID-2]];
    }
    else{
        B = vertices[faces[edgeID+1]];
    }

    return glm::length(B-A);
}

int checkForFlip(unsigned int goneVertexID, unsigned int keptVertexID, glm::vec3 newVertexPosition){
    // Loop through all the faces to find each triangle the vertex is used in
    for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
        // Find a triangle
        if(faces[edgeID] == goneVertexID){
            int triangleID = (edgeID / 3) * 3;
            glm::vec3 A = vertices[faces[triangleID]];
            glm::vec3 B = vertices[faces[triangleID + 1]];
            glm::vec3 C = vertices[faces[triangleID + 2]];
            // If the triangle will be removed skip this one
            if(faces[triangleID] == keptVertexID || faces[triangleID + 1] == keptVertexID || faces[triangleID + 2] == keptVertexID){
                continue;
            }
            else{
                // Get the normal of the triangle
                glm::vec3 oldNorm = glm::normalize(glm::cross(B-A, C-A));
                if(faces[triangleID] == goneVertexID){
                    A = newVertexPosition;
                }
                else if(faces[triangleID + 1] == goneVertexID){
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