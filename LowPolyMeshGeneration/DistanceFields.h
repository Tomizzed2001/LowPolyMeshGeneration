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

indexFaced inputMesh;

std::vector<std::vector<std::vector<glm::vec3>>> boundingBoxNodes;

std::vector<glm::vec3> fNormals;
std::vector<glm::vec3> eNormals;

float sizeOfGrid;

void fitToGrid(float *point, bool isMin);

float distanceToTriangle(glm::vec3 P, glm::vec3 A, glm::vec3 B, glm::vec3 C, unsigned int faceID);

float distanceToEdge(glm::vec3 P, glm::vec3 A, glm::vec3 B);

int sign(float value);

glm::vec3 getVertexNormal(unsigned int vertexID);

glm::vec3 findOtherHalfNormal(unsigned int v0, unsigned int v1);

//std::unordered_set<glm::vec3> getBBs(glm::vec3 point);
