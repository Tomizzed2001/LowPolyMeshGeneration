// Simplification.h
#pragma once

#include <iostream>
#include <vector>
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
std::vector<glm::mat4> quadrics;

glm::mat4 findK(int triangleID);