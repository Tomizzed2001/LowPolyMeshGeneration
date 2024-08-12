// DistanceFields.h
#pragma once

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <unordered_set>
#include <chrono>
#include "glm/glm.hpp"

struct indexFaced {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> faces;
    std::vector<glm::vec3> vNormals;
};

// Mesh data structure
indexFaced inputMesh;

// Boundding box locations
std::vector<std::vector<std::vector<glm::vec3>>> boundingBoxNodes;

// Normals
std::vector<glm::vec3> fNormals;
std::vector<glm::vec3> eNormals;

// Size of the grid
float sizeOfGrid;

/// @brief Function to round a point so that there is a nice integer number of grid points
/// @param point the point to round
/// @param isMin is it the minimum point
void fitToGrid(float *point, bool isMin);

/// @brief Function to return the distance to a triangle from a point given the inputs:
/// @param P is the point
/// @param A is a vertex on the triangle
/// @param B is a vertex on the triangle
/// @param C is a vertex on the triangle
/// @param faceID is the faceID of the triangle
/// @return distance from point to triangle
float distanceToTriangle(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C, unsigned int faceID);

/// @brief Function to return the distance to an edge from a point given the inputs:
/// @param P is the point
/// @param A is a vertex on the edge
/// @param B is a vertex on the edge
/// @return distance from point to edge
float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B);

/// @brief Get the sign of a value
/// @param value the value
/// @return the sign
int sign(float value);

/// @brief Returns the vertex normal
/// @param vertexID the vertex
/// @return the normal
glm::vec3 getVertexNormal(unsigned int vertexID);

