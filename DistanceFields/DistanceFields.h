// DistanceFields.h
#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "glm/glm.hpp"

std::vector<glm::mat4> transforms;
std::vector<std::vector<glm::vec2>> transformedFaces; // A, B, C, nAB, nBC, nCA
std::vector<glm::vec3> faceNormals;

void fitToGrid(float *point, bool isMin);

glm::mat4 getTransformMatrix(glm::vec3 A, glm::vec3 B, glm::vec3 C);

void preComputeFace(int faceID, glm::mat4 transform, glm::vec3 A, glm::vec3 B, glm::vec3 C);

float distToTriangle(int faceID, glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C);


