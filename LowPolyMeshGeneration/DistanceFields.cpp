// DistanceFields.cpp : Defines the entry point for the application.

#include "DistanceFields.h"

#include <glm/gtx/string_cast.hpp>

using namespace std;

#define GRID_SIZE 0.1

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
                faces.push_back(f);
                //cout << "Face: " << glm::to_string(faces[11]) << endl;
            }
        }
    }
    dField.close();
    
    // Calculate face normals from scratch
    fNormals.resize(faces.size());
    for(size_t fID = 0; fID < faces.size(); fID++){
        glm::vec3 AB = vertices[faces[fID].y] - vertices[faces[fID].x];
        glm::vec3 AC = vertices[faces[fID].z] - vertices[faces[fID].x];
        fNormals[fID] = glm::normalize(glm::cross(AB, AC));
    }


    // Calculate the edge normals from scratch
    eNormals.resize(faces.size()*3);
    for(size_t fID = 0; fID < faces.size(); fID++){
        // AB edge
        eNormals[fID * 3] = glm::normalize((findOtherHalfNormal(faces[fID].x, faces[fID].y) + fNormals[fID]) / float(2));
        // BC edge
        eNormals[fID * 3 + 1] = glm::normalize((findOtherHalfNormal(faces[fID].y, faces[fID].z) + fNormals[fID]) / float(2));
        // CA edge
        eNormals[fID * 3 + 2] = glm::normalize((findOtherHalfNormal(faces[fID].z, faces[fID].x) + fNormals[fID]) / float(2));
    }

    // Calculate the vertex normals from scratch
    vNormals.resize(vertices.size());
    for(size_t vID = 0; vID < vertices.size(); vID++){
        vNormals[vID] = getVertexNormal(vID);
    }

    //Round the points so the grid fits nicely
    cout << "Min Point: " << minPoint[0] << ", " << minPoint[1] << ", " << minPoint[2] << endl;
    cout << "Max Point: " << maxPoint[0] << ", " << maxPoint[1] << ", " << maxPoint[2] << endl;
    fitToGrid(minPoint, true);
    fitToGrid(maxPoint, false);
    cout << "Min Point: " << minPoint[0] << ", " << minPoint[1] << ", " << minPoint[2] << endl;
    cout << "Max Point: " << maxPoint[0] << ", " << maxPoint[1] << ", " << maxPoint[2] << endl;

    // Set up the scalar field
    int xSize = round(((maxPoint[0] - minPoint[0]) / GRID_SIZE) + 1);
    int ySize = round(((maxPoint[1] - minPoint[1]) / GRID_SIZE) + 1);
    int zSize = round(((maxPoint[2] - minPoint[2]) / GRID_SIZE) + 1);
    std::vector<std::vector<std::vector<float>>> scalarField(xSize, vector<vector<float>>(ySize, vector<float>(zSize)));

    // Get the transform matrices for each triangle
    transformedFaces.resize(faces.size()+1, vector<glm::vec2>(6));
    faceNormals.resize(faces.size()+1);
    for(size_t faceID = 0; faceID < faces.size(); faceID++){
        transforms.push_back(getTransformMatrix(vertices[faces[faceID].x], vertices[faces[faceID].y], vertices[faces[faceID].z]));
        preComputeFace(faceID, transforms.back(), vertices[faces[faceID].x], vertices[faces[faceID].y], vertices[faces[faceID].z]);
    }

    // Loop through all points in the scalar field and get the distance 
    // from them to the closest triangle
    int id;
    for (float xID = 0, x = minPoint[0]; xID < xSize; xID++, x+=GRID_SIZE) {
		for (float yID = 0, y = minPoint[1]; yID < ySize; yID++, y+=GRID_SIZE) {
			for (float zID = 0, z = minPoint[2]; zID < zSize; zID++, z+=GRID_SIZE) {
                float closestDistance = 1000;
                if( (x > -1.31 && x < -1.29) && (y > -0.51 && y < -0.49) && (z > -0.51 && z < -0.49) ){
                    cout << "Setting true" << endl;
                    here = true;
                }
                for(size_t faceID = 0; faceID < faces.size(); faceID++){
                    /*
                    float distance = distToTriangle(
                        faceID,
                        glm::vec3(x,y,z),
                        vertices[faces[faceID].x],
                        vertices[faces[faceID].y], 
                        vertices[faces[faceID].z]
                    );
                    float distance = distance3D(
                        glm::vec3(x,y,z),
                        vertices[faces[faceID].x],
                        vertices[faces[faceID].y], 
                        vertices[faces[faceID].z]
                    );
                    */
                    float distance = closestPointTriangle(
                        glm::vec3(x,y,z),
                        vertices[faces[faceID].x],
                        vertices[faces[faceID].y], 
                        vertices[faces[faceID].z],
                        faces[faceID].x,
                        faces[faceID].y,
                        faces[faceID].z,
                        faceID
                    );

                    /*
                    glm::vec3 AB = vertices[faces[faceID].y] - vertices[faces[faceID].x];
                    glm::vec3 AC = vertices[faces[faceID].z] - vertices[faces[faceID].x];
                    glm::vec3 n = glm::normalize(glm::cross(AB, AC));
                    float s = sign(glm::dot(point - glm::vec3(x,y,z), n));
                    if(s == 0){
                        s = sign(glm::dot(point - glm::vec3(x,y,z), n));
                    }
                    float distance = s * glm::length(point - glm::vec3(x,y,z));
                    */

                    if( (x > -1.31 && x < -1.29) && (y > -0.51 && y < -0.49) && (z > -0.51 && z < -0.49) ){
                        if(abs(distance) < 0.5){
                            cout << "D: " << distance << endl;
                            cout << "Vertices: " << glm::to_string(vertices[faces[faceID].x]) << glm::to_string(vertices[faces[faceID].y]) << glm::to_string(vertices[faces[faceID].z]) << endl;
                            //cout << glm::dot(point - glm::vec3(x,y,z), n) << endl;
                        }
                    }


                    if(abs(distance) < abs(closestDistance)){
                        closestDistance = distance;
                        id = faceID;
                    }
                }
                // Remove any "negative" zeros
                if(closestDistance <= 0.0001 && closestDistance >= -0.0001){
                    closestDistance = 0;
                }
                if(closestDistance > 0){
                    //cout << "New Point: " << x << ", " << y << ", " << z << " Triangle: " << id << " " << " Distance: " << closestDistance << endl;
                    //cout << "Vertices: " << glm::to_string(vertices[faces[id].x]) << glm::to_string(vertices[faces[id].y]) << glm::to_string(vertices[faces[id].z]) << endl;
                }
                scalarField[xID][yID][zID] = closestDistance;
                here = false;
            }
        }
    }

    // Dump the scalar field to a .txt file
    std::ofstream out("scalarField.txt");
    out << xSize << " " << ySize << " " << zSize << endl;
    for (int x = 0; x < xSize; x++) {
		for (int y = 0; y < ySize; y++) {
			for (int z = 0; z < zSize; z++) {
                // Truncate at 4 decimal places and add 5 as the 5th (0.00005)
				out << std::fixed << setprecision(4) << scalarField[x][y][z] << 5 << " ";
			}
            out << endl;
		}
        out << endl;
	}
    out.close();
    /*
    */


    // Testing stuff (It seems that it flips the distance so - is on the outside? Requires further testing to be
    // sure but test with real distance and can do the paper calculations too).
    //glm::mat4 T = glm::mat4(1);
    //transforms[0] = T;
    /*
    glm::vec3 A = vertices[faces[3].x];
    std::cout << glm::to_string(A) << endl;
    glm::vec3 B = vertices[faces[3].y];
    std::cout << glm::to_string(B) << endl;
    glm::vec3 C = vertices[faces[3].z];
    std::cout << glm::to_string(C) << endl;
    glm::vec3 P = glm::vec3(-2.5, 0, -2.5);
    */
    /*
    glm::vec3 A = vertices[faces[952].x];
    std::cout << glm::to_string(A) << endl;
    glm::vec3 B = vertices[faces[952].y];
    std::cout << glm::to_string(B) << endl;
    glm::vec3 C = vertices[faces[952].z];
    std::cout << glm::to_string(C) << endl;
    glm::vec3 P = glm::vec3(2, -1, -2);
    cout << "Run" << endl;
    cout << distToTriangle(3, P, A, B, C) << endl;

    glm::vec3 W = glm::vec3(1, -1, 0);
    glm::vec3 X = glm::vec3(0, 0, 0);
    glm::vec3 Y = glm::vec3(2, 0, 0);

    cout << "TIS " << dist(W,X,Y) << endl;
    */

	return 0;
}

void fitToGrid(float *point, bool isMin){
    //Get the sign
    int sign;
    if(isMin){
        sign = -2;
    }else{
        sign = 2;
    }
    for(int i = 0; i < 3; i++){
        if(point[i] < 0){
            point[i] = -(abs(point[i]) + (GRID_SIZE - fmod(abs(point[i]), GRID_SIZE))) + (sign * GRID_SIZE);
        }
        else{
            point[i] = point[i] + (GRID_SIZE - fmod(point[i], GRID_SIZE)) + (sign * GRID_SIZE);
        }
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

    glm::mat4 transform = R * T;


    // Check to see if C is in positive y or negative y
    // Flip if negative in the z axis
    if((transform * glm::vec4(C,1)).y < 0){
        glm::mat4 zRot = glm::mat4(1);
        zRot[0] = glm::vec4(-1,0,0,0);
        zRot[1] = glm::vec4(0,-1,0,0);
        R = R * zRot;
        transform = R * T;
    }


    // Check to see if B is in positive z or negative z
    // Flip if negative in the y axis
    if((transform * glm::vec4(B,1)).z < 0){
        glm::mat4 yRot = glm::mat4(1);
        yRot[0] = glm::vec4(-1,0,0,0);
        yRot[2] = glm::vec4(0,0,-1,0);
        R = R * yRot;
        transform = R * T;
    }
    /*
    */

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

    // Get the normals off each edge
    glm::vec2 AB = transformedFaces[faceID][1] - transformedFaces[faceID][0];
    transformedFaces[faceID][3] = glm::vec2(AB.y, -AB.x);
    glm::vec2 BC = transformedFaces[faceID][2] - transformedFaces[faceID][1];
    transformedFaces[faceID][4] = glm::vec2(BC.y, -BC.x);
    glm::vec2 CA = transformedFaces[faceID][0] - transformedFaces[faceID][2];
    transformedFaces[faceID][5] = glm::vec2(CA.y, -CA.x);

    // Get the normal of the face and position it off of A
    glm::vec3 n = glm::normalize(glm::cross(B-A, C-A));
    faceNormals[faceID] = A + n;
}

float distToTriangle(int faceID, glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C){
    // Calculate the normal of the triangle
    glm::vec3 ABTest = B - A;
    glm::vec3 ACTest = C - A;
    glm::vec3 normal = glm::normalize(glm::cross(ABTest, ACTest));

    // Transform P into 2D plane of triangle
    glm::vec3 newP = transforms[faceID] * glm::vec4(P,1);
    //std::cout << "newP: " << glm::to_string(newP) << endl;
    glm::vec2 tP = glm::vec2(newP.z, newP.y);
    //std::cout << "tP: " << glm::to_string(tP) << endl;
    glm::vec2 tA = transformedFaces[faceID][0];
    //std::cout << "tA: " << glm::to_string(tA) << endl;
    glm::vec2 tB = transformedFaces[faceID][1];
    //std::cout << "tB: " << glm::to_string(tB) << endl;
    glm::vec2 tC = transformedFaces[faceID][2];
    //std::cout << "tC: " << glm::to_string(tC) << endl;
    glm::vec3 n = transforms[faceID] * glm::vec4(faceNormals[faceID],1);
    //cout << "Norm: " << glm::to_string(n) << endl;

    // Determine if point is within the bounds of the triangle using the half plane test
    float AB = (tB.y - tA.y) * tP.x + (tA.x - tB.x) * tP.y + (tB.x * tA.y - tA.x * tB.y);
    //std::cout << "AB: " << AB << endl;
    float BC = (tC.y - tB.y) * tP.x + (tB.x - tC.x) * tP.y + (tC.x * tB.y - tB.x * tC.y);
    //std::cout << "BC: " << BC << endl;
    float CA = (tA.y - tC.y) * tP.x + (tC.x - tA.x) * tP.y + (tA.x * tC.y - tC.x * tA.y);
    //std::cout << "CA: " << CA << endl;
    if(AB <= 0 && BC <= 0 && CA <= 0){
        //std::cout << "INSIDE" << endl;
        // Negative when behind the triangle
        return n.x * newP.x;
    }

    // Point must be outside the triangle so check to see if there is a closest vertex
    glm::vec2 AnAB = tA + transformedFaces[faceID][3];
    glm::vec2 AnCA = tA + transformedFaces[faceID][5];
    glm::vec2 BnAB = tB + transformedFaces[faceID][3];
    glm::vec2 BnBC = tB + transformedFaces[faceID][4];
    glm::vec2 CnBC = tC + transformedFaces[faceID][4];
    glm::vec2 CnCA = tC + transformedFaces[faceID][5];

    // Test for A vertex
    float dAnAB = (AnAB.y - tA.y) * tP.x + (tA.x - AnAB.x) * tP.y + (AnAB.x * tA.y - tA.x * AnAB.y);
    float dAnCA = (AnCA.y - tA.y) * tP.x + (tA.x - AnCA.x) * tP.y + (AnCA.x * tA.y - tA.x * AnCA.y);
    if(dAnAB >= 0 && dAnCA <= 0){
        //cout << "A" << endl;
        float s = sign(glm::dot(P-A, normal));
        return s * glm::length(P-A);
    }

    // Test for B vertex
    float dBnAB = (BnAB.y - tB.y) * tP.x + (tB.x - BnAB.x) * tP.y + (BnAB.x * tB.y - tB.x * BnAB.y);
    float dBnBC = (BnBC.y - tB.y) * tP.x + (tB.x - BnBC.x) * tP.y + (BnBC.x * tB.y - tB.x * BnBC.y);
    if(dBnAB <= 0 && dBnBC >= 0){
        //cout << "B" << endl;
        float s = sign(glm::dot(P-B, normal));
        return s * glm::length(P-B);
    }

    // Test for C vertex
    float dCnBC = (CnBC.y - tC.y) * tP.x + (tC.x - CnBC.x) * tP.y + (CnBC.x * tC.y - tC.x * CnBC.y);
    float dCnCA = (CnCA.y - tC.y) * tP.x + (tC.x - CnCA.x) * tP.y + (CnCA.x * tC.y - tC.x * CnCA.y);
    if(dCnBC <= 0 && dCnCA >= 0){
        //cout << "C" << endl;
        float s = sign(glm::dot(P-C, normal));
        return s * glm::length(P-C);
    }

    // Test for AB
    if(AB > 0){
        //std::cout << "AB" << endl;
        return distanceToEdge(P, A, B, normal);
    }
    // Test for BC
    if(BC > 0){
        //cout << "BC" << endl;
        return distanceToEdge(P, B, C, normal);
    }
    // Test for CA
    if(CA > 0){
        //cout << "CA" << endl;
        return distanceToEdge(P, C, A, normal);
    }

    return 0;
    
}

float distance3D(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C){
    // Check that the point does not lie on one of the vertices
    if(P == A || P == B || P == C){
        return 0;
    }

    // Calculate the normal of the triangle
    glm::vec3 AB = B - A;
    glm::vec3 AC = C - A;
    glm::vec3 n = glm::normalize(glm::cross(AB, AC));

    // Angle between normal and AP
    glm::vec3 AP = glm::normalize(P - A);
    float alpha = glm::dot(AP, n);

    // Vector P to P projected on the plane
    glm::vec3 PA = A - P;
    float lengthPP = glm::length(PA) * alpha;
    glm::vec3 PP = -lengthPP * n;

    // Projection of P is then
    glm::vec3 projP = P + PP;
    cout << glm::to_string(projP) << endl;

    // Determine if it is inside the triangle
    


    // Can be pre-calculated
    glm::vec3 V1 = glm::normalize(A-B) + glm::normalize(A-C);
    cout << "V1: " << glm::to_string(V1) << endl;
    glm::vec3 V2 = glm::normalize(B-C) + glm::normalize(B-A);
    cout << "V2: " << glm::to_string(V2) << endl;
    glm::vec3 V3 = glm::normalize(C-A) + glm::normalize(C-B);
    cout << "V3: " << glm::to_string(V3) << endl;

    float f1 = glm::dot(glm::cross(V1, projP-A), n);
    float f2 = glm::dot(glm::cross(V2, projP-B), n);
    float f3 = glm::dot(glm::cross(V3, projP-C), n);





    // Determine if point is within the bounds of the triangle using the half plane test
    float testAB = (B.y - A.y) * projP.x + (A.x - B.x) * projP.y + (B.x * A.y - A.x * B.y);
    cout << testAB << endl;
    float testBC = (C.y - B.y) * projP.x + (B.x - C.x) * projP.y + (C.x * B.y - B.x * C.y);
    cout << testBC << endl;
    float testCA = (A.y - C.y) * projP.x + (C.x - A.x) * projP.y + (A.x * C.y - C.x * A.y);
    cout << testCA << endl;
    if(testAB >= 0 && testBC >= 0 && testCA >= 0){
        cout << "INSIDE" << endl;
        float s = sign(glm::dot(P-projP, n));
        return s * lengthPP;
    }

    // Must be outside the triangle
    // Closest point is A
    if(testAB <=0 && testCA <=0){
        cout << "A" << endl;
        float s = sign(glm::dot(P-A, n));
        return s * glm::length(P-A);
    }
    // Closest point is B
    if(testAB <=0 && testBC <=0){
        cout << "B" << endl;
        float s = sign(glm::dot(P-B, n));
        return s * glm::length(P-B);
    }
    // Closest point is C
    if(testBC <=0 && testCA <=0){
        cout << "C" << endl;
        float s = sign(glm::dot(P-C, n));
        return s * glm::length(P-C);
    }

    // Closest point lies on AB
    if(testAB <= 0){
        cout << "AB" << endl;
        return distanceToEdge(P, A, B, n);
    }
    // Closest point lies on BC
    if(testBC <= 0){
        cout << "BC" << endl;
        return distanceToEdge(P, B, C, n);
    }
    // Closest point lies on CA
    if(testCA <= 0){
        cout << "CA" << endl;
        return distanceToEdge(P, C, A, n);
    }

    // Return 0 to prevent error warnings
    return 0;
}

float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 N){
    glm::vec3 AB = B - A;
    // Project P onto the line
    float projection = clampit(((P.x - A.x) * (B.x - A.x) + (P.y - A.y) * (B.y - A.y) + (P.z - A.z) * (B.z - A.z)) / (AB.x*AB.x + AB.y*AB.y + AB.z*AB.z), 0, 1);
    glm::vec3 projectedPoint = glm::vec3(A.x + projection * (B.x - A.x), A.y + projection * (B.y - A.y), A.z + projection * (B.z - A.z));
    float s = sign(glm::dot(P-projectedPoint, N));
    return s * glm::length(P - projectedPoint);
}

float dist(glm::vec3 P, glm::vec3 A, glm::vec3 B){
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

float dot2(glm::vec3 v){
    return glm::dot(v,v);
}

float clampit(float v, float v1, float v2){
    if(v < v1){
        return v1;
    }
    else if(v > v2){
        return v2;
    }
    else{
        return v;
    }
}

float closestPointTriangle(glm::vec3 const& p, glm::vec3 const& a, glm::vec3 const& b, glm::vec3 const& c, unsigned int aID, unsigned int bID, unsigned int cID, unsigned int fID)
{
    int mode = 0;

    const glm::vec3 ab = b - a;
    const glm::vec3 ac = c - a;
    const glm::vec3 ap = p - a;

    const float d1 = glm::dot(ab, ap);
    const float d2 = glm::dot(ac, ap);
    // If vertex A
    if (d1 <= 0.f && d2 <= 0.f){
        mode = 1;
        if(here){
            cout << mode << endl;
        }

        glm::vec3 point = a;

        // Calculate the sign using vertex normal
        glm::vec3 n = vNormals[aID];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    } 
        

    const glm::vec3 bp = p - b;
    const float d3 = glm::dot(ab, bp);
    const float d4 = glm::dot(ac, bp);
    //If vertex B
    if (d3 >= 0.f && d4 <= d3){
        mode = 2;
        if(here){
            cout << mode << endl;
        }

        glm::vec3 point = b; //#2

        // Calculate the sign using vertex normal
        glm::vec3 n = vNormals[bID];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    } 

    const glm::vec3 cp = p - c;
    const float d5 = glm::dot(ab, cp);
    const float d6 = glm::dot(ac, cp);
    // If vertex C
    if (d6 >= 0.f && d5 <= d6){
        mode = 3;
        if(here){
            cout << mode << endl;
        }

        glm::vec3 point = c; //#3

        // Calculate the sign using vertex normal
        glm::vec3 n = vNormals[cID];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    } 

    const float vc = d1 * d4 - d3 * d2;
    // If edge AB
    if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
    {
        // DEBUG LINES
        mode = 4;
        if(here){
            cout << mode << endl;
        }

        const float v = d1 / (d1 - d3);
        glm::vec3 point = a + v * ab;
        
        // Calculate the sign using triangle normal
        glm::vec3 n = eNormals[fID * 3];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    }
        
    const float vb = d5 * d2 - d1 * d6;
    // If edge CA
    if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
    {
        mode = 5;
        if(here){
            cout << mode << endl;
        }

        const float v = d2 / (d2 - d6);
        glm::vec3 point = a + v * ac; //#5

        // Calculate the sign using triangle normal
        //glm::vec3 n = fNormals[fID];
        glm::vec3 n = eNormals[fID * 3 + 2];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    }
        
    const float va = d3 * d6 - d5 * d4;
    // If edge BC
    if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
    {
        mode = 6;
        if(here){
            cout << mode << endl;
        }

        const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        glm::vec3 point = b + v * (c - b); //#6

        // Calculate the sign using triangle normal
        //glm::vec3 n = fNormals[fID];
        glm::vec3 n = eNormals[fID * 3 + 1];
        float s = sign(glm::dot(point - p, n));
        float d = glm::length(point - p);

        return s * d;
    }

    // Must be inside
    mode = 0;
    if(here){
        cout << mode << endl;
    }

    const float denom = 1.f / (va + vb + vc);
    const float v = vb * denom;
    const float w = vc * denom;

    glm::vec3 point = a + v * ab + w * ac; //#0

    // Calculate the sign using triangle normal
    //glm::vec3 n = fNormals[fID];
    glm::vec3 n = fNormals[fID];
    float s = sign(glm::dot(point - p, n));
    float d = glm::length(point - p);

    return s * d;
}

glm::vec3 getVertexNormal(unsigned int vertexID){
    vector<glm::vec3> triangleNormals;
    // Find the triangles using the vertex
    for(size_t fID = 0; fID < faces.size(); fID++){
        if(faces[fID].x == vertexID || faces[fID].y == vertexID || faces[fID].z == vertexID){
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
    for(size_t fID = 0; fID < faces.size(); fID++){
        if((faces[fID].x == v1 && faces[fID].y == v0) || 
            (faces[fID].y == v1 && faces[fID].z == v0) ||
            (faces[fID].z == v1 && faces[fID].x == v0)
        ){
            return fNormals[fID];
        }
    }
    return glm::vec3(0);
}