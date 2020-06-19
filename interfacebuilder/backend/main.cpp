#include "Python.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "logging.h"
#include "backend.h"

using float_2d_vec = std::vector<std::vector<float>>;

int main()
{
    std::vector<float> a1{3.5, 0};
    std::vector<float> a2{0, 7};
    std::vector<float> b1{3.5, 0};
    std::vector<float> b2{0, 7};
    std::vector<float> angles{0, 1, 2, 3, 4, 5};
    int Ntrans = 3;
    float crit = 0.05;
    float mingamma = 20;
    float maxgamma = 160;
    float_2d_vec pairs;
    pairs = backend_routine(a1, a2, b1, b2, angles, Ntrans, crit, mingamma, maxgamma);
    log_2dvec(pairs);
}
