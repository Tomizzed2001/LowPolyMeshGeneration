// DistanceFields.h
#pragma once

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
#include "glm/glm.hpp"

std::vector<glm::vec3> vertices; // x y z
std::vector<glm::vec3> faces;
std::vector<glm::vec3> vNormals;
std::vector<glm::vec3> fNormals;
std::vector<glm::vec3> eNormals;

std::vector<glm::mat4> transforms;
std::vector<std::vector<glm::vec2>> transformedFaces; // A, B, C, nAB, nBC, nCA
std::vector<glm::vec3> faceNormals;
std::vector<glm::vec3> newFaceNormals;

bool here = false;

void fitToGrid(float *point, bool isMin);

float distanceToTriangle(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C, unsigned int faceID);

float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B);

int sign(float value);

glm::vec3 getVertexNormal(unsigned int vertexID);

glm::vec3 findOtherHalfNormal(unsigned int v0, unsigned int v1);
