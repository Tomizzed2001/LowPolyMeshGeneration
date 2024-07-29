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
#include <chrono>
#include "glm/glm.hpp"

struct dirEdge {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> vertexNormals;
    std::vector<unsigned int> edges;
    std::vector<int> otherhalves;
};

dirEdge mesh;

std::vector<std::unordered_set<unsigned int>> oneRings;

std::vector<unsigned int> firstDirectedEdges;

std::vector<glm::mat4> quadrics;    // Q matrix per vertex
std::vector<float> errorCosts;      // Cost for each half edge
std::vector<glm::vec3> optimalVertexPosition;   // Per vertex

std::unordered_set<int> updatedEdges;

glm::mat4 findK(int triangleID);

std::unordered_set<unsigned int> findOneRing(unsigned int vertexID);

void removeFace(unsigned int faceID, int faceNum);

void findOtherHalf(unsigned int edgeID);

void updateQ(unsigned int vertexID);

float getEdgeError(unsigned int edgeID);

float getEdgeLength(unsigned int edgeID);

void outputToDiredge();

void outputToObject();

void outputToObject(int num);

int checkForFlip(unsigned int goneVertexID, unsigned int keptVertexID, glm::vec3 newVertexPosition);