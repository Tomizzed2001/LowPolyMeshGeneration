// MarchingCubes.h
#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "glm/glm.hpp"

std::vector<std::vector<std::vector<float>>> scalarField;

std::vector<glm::vec3> vertices;
std::vector<int> faces;
std::vector<int> firstDirectedEdges;
std::vector<int> otherHalf;

int getVertexID(glm::vec3 vertex);