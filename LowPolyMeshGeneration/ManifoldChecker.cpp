#include "ManifoldChecker.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

int main(int argc, char** argv){
    //Check the inputs
    //Check there is exactly one input file and return error if not
    if (argc != 2){
        std::cout << "There is " << argc -1 << " input files program needs exactly one .face file" << std::endl;
        return 0;
    }

    //Check the file is a .diredge file
    
    std::string s = argv[1];
    if (s.find(".diredge") != string::npos){
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
                    mesh.firstDirectedEdges.push_back(id);
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
    }
    else if(s.find(".obj") != string::npos){
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
                    }
                    inputMesh.vertices.push_back(p);
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
                        // Ignore the normal
                        getline(ssFace, vtn, '/');
                    }
                    inputMesh.faces.push_back(f);
                    //cout << "Face: " << glm::to_string(faces[11]) << endl;
                }
            }
        }
        dField.close();

        mesh.vertices = inputMesh.vertices;
        for(size_t i = 0; i < inputMesh.faces.size(); i++){
            mesh.edges.push_back(int(inputMesh.faces[i][0]));
            mesh.edges.push_back(int(inputMesh.faces[i][1]));
            mesh.edges.push_back(int(inputMesh.faces[i][2]));
        }

        // Index storage variables
        unsigned int v0, v1;

        mesh.otherhalves.resize(mesh.edges.size(),-1);

        // Get the other half
        for(size_t eID = 0; eID < mesh.edges.size(); eID++){
            // Since faces form the edges check if it needs to search for 2 to 0 edge
            if (eID % 3 == 2){
                v0 = mesh.edges[eID];
                v1 = mesh.edges[eID-2];
            }
            else{
                v0 = mesh.edges[eID];
                v1 = mesh.edges[eID+1];
            }

            // Look for v1 to v0
            for(size_t fID = 0; fID < mesh.edges.size(); fID+=3){
                if(mesh.edges[fID] == v1 && mesh.edges[fID+1] == v0){
                    mesh.otherhalves[eID] = fID + 0;
                    break;
                }
                else if(mesh.edges[fID+1] == v1 && mesh.edges[fID+2] == v0){
                    mesh.otherhalves[eID] = fID + 1;
                    break;
                }
                else if(mesh.edges[fID+2] == v1 && mesh.edges[fID] == v0){
                    mesh.otherhalves[eID] = fID + 2;
                    break;
                }
            }
        }
    }
    else{
        std::cout << "Program only takes input of a .diredge or .obj file" << std::endl;
        return 0;
    }

    //Test for manifold
    bool isManifold = true;

    //Check that each edge occurs exactly twice otherwise hole
    for(size_t edgeID = 0; edgeID < mesh.otherhalves.size(); edgeID++){
        if(int(edgeID) != mesh.otherhalves[mesh.otherhalves[edgeID]]){
            cout << "Not 2-Manifold. Edge: " << edgeID << " does not share exactly 2 faces." << endl;
            isManifold = false;
        }
    }
    if(!isManifold) return 0;

    //Make the dual graph (nodes are faces, edges are shared edges between the faces)
    //Assign empty variables and arrays
    neighbourTriangles.resize(mesh.edges.size() / 3, vector<int>(3, int()));
    //For each face find all neighbour faces
    for(size_t faceID = 0; faceID < (mesh.edges.size() / 3); faceID++){
        for(size_t i = 0; i < 3; i++){
            neighbourTriangles[faceID][i] = mesh.otherhalves[faceID*3+i] / 3;
        }
    }

    //Find all faces that can be reached from face 0 using BFS
    BFS(0);

    int faceNum = (mesh.edges.size() / 3);

    //Check to see if all faces are present in the dual graph
    int unfoundFace, facesFound = 0;

    for(int i = 0; i < faceNum; i++){
        //If the value is 0 the face was never found
        if(dualGraph[i] != 1){
            //The face has not been found so must be in another object or beyond pinch point
            unfoundFace = i;
        }
        //Face 
        else{
            facesFound++;
        }
    }

    //This is only called if there is more than one object
    //This section find the vertex of the pinch point
    if(facesFound != faceNum){
        //Max number of objects is number of faces / 4 since smallest manifold object needs 4 faces
        vector<vector<int>> objectVertices;
        objectVertices.resize(faceNum/4, vector<int>(mesh.vertices.size(), int()));

        //Array for storing the faces visited
        vector<int> visitedFaces;
        visitedFaces.resize(faceNum);
        //Currently working on the first object
        int objectNum = 0;
        //Put the dual graphs vertices in the place of the corresponding object
        for(int dualID = 0; dualID < faceNum; dualID++){
            //Face is present
            if (dualGraph[dualID] == 1){
                visitedFaces[dualID] = 1;
                objectVertices[0][ mesh.edges[dualID * 3]     ] = 1;
                objectVertices[0][ mesh.edges[dualID * 3 + 1] ] = 1;
                objectVertices[0][ mesh.edges[dualID * 3 + 2] ] = 1;
            }
        }

        //Find all other objects in the mesh
        while(facesFound != faceNum){
            //Next object
            objectNum++;
            //Do BFS to find another set of faces
            BFS(unfoundFace);
            
            //Put the dual graphs vertices in the place of the corresponding object
            for(int dualID = 0; dualID < faceNum; dualID++){
                //Face is present
                if (dualGraph[dualID] == 1){
                    visitedFaces[dualID] = 1;
                    objectVertices[objectNum][ mesh.edges[dualID * 3]     ] = 1;
                    objectVertices[objectNum][ mesh.edges[dualID * 3 + 1] ] = 1;
                    objectVertices[objectNum][ mesh.edges[dualID * 3 + 2] ] = 1;
                }
            }
            
            //Increment the number of new faces found 
            for(int i = 0; i < faceNum; i++){
                //If the value is 0 the face was never found
                if(dualGraph[i] != 1 && visitedFaces[i] == 0){
                    //The face has not been found so must be in another object or beyond pinch point
                    unfoundFace = i;
                }
                else if(dualGraph[i] == 1){
                    facesFound++;
                }
            }
        }

        vector<int> commonVertices(mesh.vertices.size());
        //Compare the vertices of each object to check for a common vertex (pinch point)
        for(size_t vID = 0; vID < mesh.vertices.size(); vID++){
            for(int objID = 0; objID <= objectNum; objID++){
                commonVertices[vID] += objectVertices[objID][vID];
                if(commonVertices[vID] > 1){
                    std::cout << "Pinch point found at vertex: " << vID << std::endl;
                    isManifold = false;
                    return 0;
                }
            }
        }
        
        //If there is multiple objects and no pinch point then it is not manifold
        if (objectNum > 0 && isManifold){
            std::cout << "Mesh has " << objectNum+1 << " objects, so it is not a polyhedron. Therefore, it cannot be manifold. " << std::endl;
            isManifold = false;
            return 0;
        }
    }



    if(isManifold){
        cout << "Is 2-Manifold!" << endl;
    }

}

void BFS(int startingNode)
{
 
    //Initialise arrays and queue indexes
    dualGraph.resize(mesh.edges.size()/3);
    vector<int> queue;
    queue.resize(mesh.edges.size()/3);
    int currentIndex = 0, endIndex = 0, currentFace = startingNode;

    //Start with face 0 marked as seen and add to the queue
    dualGraph[startingNode] = 1;
    queue[endIndex++] = startingNode;

    while (currentIndex != endIndex){
        //Set current face and move front of the queue up one
        currentFace = queue[currentIndex++];

        //Loop through each neighbour face
        for(int i = 0; i < 3; i++){
            
            //Check if it has already been visited, if not add to queue and mark as visited
            if (dualGraph[neighbourTriangles[currentFace][i]] == 0){
                queue[endIndex++] = neighbourTriangles[currentFace][i];
                dualGraph[neighbourTriangles[currentFace][i]] = 1;

            }
        }
    }
    
    return;
}