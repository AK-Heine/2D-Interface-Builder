#pragma once
#include <vector>

using float_2d_vec = std::vector<std::vector<float>>;

float_2d_vec find_coincidence(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, std::vector<float> angles, int Ntrans, float crit);
float_2d_vec backend_routine(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, std::vector<float> angles, int Ntrans, float crit, float mingamma, float maxgamma);
