#include "Simplification.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define DESIRED_TRIANGLE_COUNT 40

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

    //return 0;

    // Calculate the Quadric Error Metric for each half-edge
    // Step 1: Find Q for each vertex
    quadrics.resize(vertices.size());
    for(size_t vID = 0; vID < vertices.size(); vID++){
        updateQ(vID);
    }
    
    // Step 2: Find error cost for each valid half edge
    errorCosts.resize(otherHalf.size());
    for(size_t ohID = 0; ohID < otherHalf.size(); ohID++){
        updateError(ohID);
    }

    // Step 3: Place otherhalf IDs on a heap ordererd by error cost.
    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> collapseOrder;
    for(size_t i = 0; i < errorCosts.size(); i++){
        if(errorCosts[i] == -1){
            continue;
        }
        collapseOrder.push(make_pair(errorCosts[i], i));
    }

    priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>> testOrder = collapseOrder;

    // DEBUG COLLAPSE ORDER / ERROR
    ofstream out("errors/error0.error");
    while(! testOrder.empty()){
        out << testOrder.top().first << " " << testOrder.top().second << endl;
        testOrder.pop();
    }
    out.close();
    // DEBUG OTHERHALVES
    std::string var = "otherHalves/otherHalf0.oh";
    ofstream outt(var);
    for(size_t i = 0; i < otherHalf.size(); i++){
        outt << otherHalf[i] << endl;
    }
    outt.close();

    vector<unsigned int> removedVertices;

    int counter = 1;

    // Complete the collapse in the order of lowest error cost first
    while((faces.size() / 3) > DESIRED_TRIANGLE_COUNT || collapseOrder.empty()){
        cout << "Run: " << counter << endl;
        // Get the edge to collapse
        unsigned int edge = collapseOrder.top().second;
        unsigned int otherEdge = otherHalf[collapseOrder.top().second];

        //cout << "Edge: " << edge << " Other Edge: " << otherEdge << endl;

        // Get the vertices involved
        unsigned int keptVertex = faces[edge]; 
        unsigned int goneVertex = faces[otherEdge];

        //cout << "Keep Vertex: " << keptVertex + 1 << " Remove Vertex: " << goneVertex + 1 << endl;

        // Check that the second triangle to be removed is not the final face
        // Remove the faces
        if ((otherEdge / 3) * 3 == faces.size()-3){
            //cout << "OTHER IS ON THE FINAL FACE" << endl;
            removeFace(otherEdge / 3);
            removeFace(edge / 3);
        }
        else{
            removeFace(edge / 3);
            removeFace(otherEdge / 3);
        }

        // Replace all instances of the goneVertex with the keptVertex
        for(size_t eID = 0; eID < faces.size(); eID++){
            if(faces[eID] == goneVertex){
                faces[eID] = keptVertex;
            }
        }

        removedVertices.push_back(goneVertex);


        // Update half edges to be merged edges
        // Final 6 half edges no longer exist
        unsigned int newEdgeNum = otherHalf.size() - 6;
        otherHalf.erase(otherHalf.end()-6, otherHalf.end());
        for(unsigned int eID = 0; eID < newEdgeNum; eID++){
            /*
            if(otherHalf[eID] >= newEdgeNum){
                findOtherHalf(eID);
            }
            */
            findOtherHalf(eID);
        }

        // Update the vertices in the new 1-ring of keptVertex
        std::unordered_set<unsigned int> aOneRing = findOneRing(keptVertex);
        updateQ(keptVertex);
        for(auto id: aOneRing) {
            updateQ(id);
        }

        // Update the error costs for any edge involving keptVertex
        vector<unsigned int> changedEdges;
        for(unsigned int eID = 0; eID < faces.size(); eID++){
            if(faces[eID] != keptVertex){
                continue;
            }
            else{
                if (eID % 3 == 0){
                    changedEdges.push_back(eID);
                    changedEdges.push_back(eID + 2);
                }
                else{
                    changedEdges.push_back(eID);
                    changedEdges.push_back(eID-1);
                }
            }
        }
        for(unsigned int eID = 0; eID < changedEdges.size(); eID++){
            updateError(changedEdges[eID]);
        }

        // Re-Generate Collapse Order
        collapseOrder = priority_queue<pair<float, int>, vector<pair<float, int>>, greater<pair<float, int>>>();
        for(size_t i = 0; i < errorCosts.size(); i++){
            if(errorCosts[i] == -1){
                continue;
            }
            collapseOrder.push(make_pair(errorCosts[i], i));
        }

        // DEBUG OTHERHALVES
        std::string var = "otherHalves/otherHalf" + to_string(counter) + ".oh";
        ofstream outt(var);
        for(size_t i = 0; i < otherHalf.size(); i++){
            outt << otherHalf[i] << endl;
        }
        outt.close();
        // DEBUG COLLAPSE ORDER / ERROR
        testOrder = collapseOrder;
        var = "errors/error" + to_string(counter) + ".error";
        ofstream out(var);
        while(! testOrder.empty()){
            out << testOrder.top().first << " " << testOrder.top().second << endl;
            testOrder.pop();
        }
        out.close();
        // DEBUG OBJECT FILE
        outputToObject(counter);
        counter++;

        /*
        if(counter == 18){
            break;
        }
        */
    }

    //outputToObject();
    //cout << "OUTPUT" << endl;

    //cout << "Removing vertices" << endl;
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

    // Output to .obj file
    outputToObject();
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
            if(edgeID % 3 == 0){
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
        vector<pair<unsigned int, unsigned int>> swapMap;
        // Swap vertices
        swapMap.push_back(make_pair(faces[faceID], faces[faces.size()-3]));
        faces[faceID] = faces[faces.size()-3];
        swapMap.push_back(make_pair(faces[faceID+1], faces[faces.size()-2]));
        faces[faceID+1] = faces[faces.size()-2];
        swapMap.push_back(make_pair(faces[faceID+2], faces[faces.size()-1]));
        faces[faceID+2] = faces[faces.size()-1];

        // Update otherHalf to match swapped faces
        for(pair<unsigned int, unsigned int> swap: swapMap){
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
    vector<glm::mat4> allK;
    glm::mat4 Q = glm::mat4(0);

    // Find the 1-ring
    for(size_t edgeID = 0; edgeID < faces.size(); edgeID++){
        if(faces[edgeID] == vertexID){
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

void updateError(unsigned int edgeID){
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
        errorCosts[edgeID] = -1;
        return;
    }
    
    // Place on the "from" vertex
    glm::vec4 v1 = glm::vec4(vertices[faces[edgeID]], 1);

    // Get the quadric for each vertex
    glm::mat4 Q1 = quadrics[faces[edgeID]];
    glm::mat4 Q2 = quadrics[faces[otherHalf[edgeID]]];

    // Join the quadrics and calculate the error vQv
    errorCosts[edgeID] = glm::length( v1 * (Q1 + Q2) * v1 );
    //glm::vec4 quadraticError =  v1 * (Q1 + Q2) * v1;
    //errorCosts[edgeID] = (quadraticError.x * v1.x) + (quadraticError.y * v1.y) + (quadraticError.z * v1.z) + quadraticError[3];

    return;
}

void outputToDiredge(){
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