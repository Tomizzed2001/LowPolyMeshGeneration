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

void fitToGrid(float *point, bool isMin);

glm::mat4 getTransformMatrix(glm::vec3 A, glm::vec3 B, glm::vec3 C);


