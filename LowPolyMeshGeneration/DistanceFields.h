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

std::vector<glm::mat4> transforms;
std::vector<std::vector<glm::vec2>> transformedFaces; // A, B, C, nAB, nBC, nCA
std::vector<glm::vec3> faceNormals;
std::vector<glm::vec3> newFaceNormals;

bool here = false;

void fitToGrid(float *point, bool isMin);

glm::mat4 getTransformMatrix(glm::vec3 A, glm::vec3 B, glm::vec3 C);

void preComputeFace(int faceID, glm::mat4 transform, glm::vec3 A, glm::vec3 B, glm::vec3 C);

float distToTriangle(int faceID, glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C);

float distance3D(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C);

float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 N);

float dist(glm::vec3 P, glm::vec3 A, glm::vec3 B);

int sign(float value);

float clampit(float v, float v1, float v2);

glm::vec3 closestPointTriangle(glm::vec3 const& p, glm::vec3 const& a, glm::vec3 const& b, glm::vec3 const& c);
