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
#include "glm/glm.hpp"

std::vector<glm::vec3> vertices;
std::vector<unsigned int> faces;
std::vector<unsigned int> firstDirectedEdges;
std::vector<unsigned int> otherHalf;
std::vector<glm::mat4> quadrics;    // Q matrix per vertex
std::vector<float> errorCosts;

glm::mat4 findK(int triangleID);

std::unordered_set<unsigned int> findOneRing(unsigned int vertexID);

void removeFace(unsigned int faceID);

void findOtherHalf(unsigned int edgeID);
