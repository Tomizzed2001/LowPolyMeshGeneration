// ManifoldChecker.h
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

struct dirEdge {
    std::vector<glm::vec3> vertices;
    std::vector<unsigned int> edges;
    std::vector<unsigned int> firstDirectedEdges;
    std::vector<int> otherhalves;
};

struct indexFaced {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> faces;
};

// Mesh data structures
dirEdge mesh;
indexFaced inputMesh;

std::vector<int> dualGraph;
std::vector<std::vector<int>> neighbourTriangles;

void BFS(int startingNode);