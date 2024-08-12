// MarchingCubes.h
#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include "glm/glm.hpp"

// Lookup table for the edge table
int cubeEdgeLookUp[12][4] = {
	{ 0,  0,  0, 1},    // 0
	{ 0,  1,  0, 2},    // 1
    { 0,  0,  1, 1},    // 2
    { 0,  0,  0, 2},    // 3
    { 1,  0,  0, 1},    // 4
    { 1,  1,  0, 2},    // 5
    { 1,  0,  1, 1},    // 6
    { 1,  0,  0, 2},    // 7
    { 0,  0,  0, 0},    // 8
    { 0,  1,  0, 0},    // 9
    { 0,  1,  1, 0},    // 10
    { 0,  0,  1, 0}     // 11
    };

// Directed edge data structure
struct dirEdge {
    std::vector<glm::vec3> vertices;
    std::vector<int> edges;
    std::vector<int> otherhalves;
};

/// Mesh data structure
dirEdge mesh;

// Input distance field
std::vector<std::vector<std::vector<float>>> scalarField;

// Edges of the entire grid
std::vector<std::vector<std::vector<std::array<int,3>>>> cubeEdges;

// Triangle soup
std::vector<glm::vec3> triSoup;

/// @brief Check if a vertex already exists and return the id of it if so
/// @param vertex The vertex position
/// @return The vertex ID
int getVertexID(glm::vec3 vertex);

/// @brief Output to a .OBJ file
void outputToObject();