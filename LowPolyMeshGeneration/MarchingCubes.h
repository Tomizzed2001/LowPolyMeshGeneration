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

std::vector<std::vector<std::vector<float>>> scalarField;

std::vector<glm::vec3> vertices;
std::vector<int> faces;
std::vector<int> firstDirectedEdges;
std::vector<int> otherHalf;

std::vector<std::vector<std::vector<std::array<int,3>>>> cubeEdges;

std::vector<glm::vec3> triSoup;

int getVertexID(glm::vec3 vertex);

void outputToObject(int num);

void outputToObject();

void getSingleCube(int caseNum);