// Simplification.h
#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <chrono>
#include "glm/glm.hpp"

// Directed edge data structure
struct dirEdge {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> vertexNormals;
    std::vector<unsigned int> edges;
    std::vector<int> otherhalves;
};

// Mesh data structue
dirEdge mesh;

// Set of vertex one rings
std::vector<std::unordered_set<unsigned int>> oneRings;

// First directed edge
std::vector<unsigned int> firstDirectedEdges;

std::vector<glm::mat4> quadrics;    // Q matrix per vertex
std::vector<float> errorCosts;      // Cost for each half edge
std::vector<glm::vec3> optimalVertexPosition;   // Per vertex

// Edges to update
std::unordered_set<int> updatedEdges;

/// @brief Return the fundamental error quadric K
/// @param triangleID The ID of the triangle
/// @return K
glm::mat4 findK(int triangleID);

/// @brief Find the one ring of a given vertex
/// @param vertexID The vertex
/// @return The one ring
std::unordered_set<unsigned int> findOneRing(unsigned int vertexID);

/// @brief Removes a face from the data structure
/// @param faceID The ID of the face to remove
/// @param faceNum The first (1) or second (2) face to be removed
void removeFace(unsigned int faceID, int faceNum);

/// @brief Find the Q value for a vertex
/// @param vertexID The vertex ID
void updateQ(unsigned int vertexID);

/// @brief Get the error of collapsing a certain edge
/// @param edgeID The edge ID
/// @return The error cost
float getEdgeError(unsigned int edgeID);

/// @brief Output to a directed edge file
void outputToDiredge();

/// @brief Output to a .OBJ file
void outputToObject();

/// @brief Check to see if a collapse will cause a triangle flip
/// @param goneVertexID Vertex to remove
/// @param keptVertexID Vertex to keep
/// @param newVertexPosition Position of the new vertex
/// @return returns -1 if invalid 0 otherwise
int checkForFlip(unsigned int goneVertexID, unsigned int keptVertexID, glm::vec3 newVertexPosition);