#include "Python.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <cmath>
#include <vector>

namespace py = pybind11;

float coincidence(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, int &m1, int &m2, int &n1, int &n2, float &theta)
{
    float v11 = a1[0] * m1 + a2[0] * m2;
    float v12 = a1[1] * m1 + a2[1] * m2;
    float v21 = std::cos(theta) * (b1[0] * n1 + b2[0] * n2) - std::sin(theta) * (b1[1] * n1 + b2[1] * n2);
    float v22 = std::sin(theta) * (b1[0] * n1 + b2[0] * n2) + std::cos(theta) * (b1[1] * n1 + b2[1] * n2);
    float norm = std::sqrt((v11 - v21) * (v11 - v21) + (v12 - v22) * (v12 - v22));
    return norm;
};

std::vector<std::vector<float>> find_coincidence(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, std::vector<float> angles, int Ntrans, float crit)
{
    std::vector<std::vector<float>> diffs;
    for (int m1 = -Ntrans; m1 < Ntrans + 1; m1++)
    {
        for (int m2 = -Ntrans; m2 < Ntrans + 1; m2++)
        {
            for (int n1 = -Ntrans; n1 < Ntrans + 1; n1++)
            {
                for (int n2 = -Ntrans; n2 < Ntrans + 1; n2++)
                {
                    for (float theta : angles)
                    {
                        float d = coincidence(a1, a2, b1, b2, m1, m2, n1, n2, theta);
                        if (d < crit)
                        {
                            std::vector<float> diff{static_cast<float>(m1), static_cast<float>(m2), static_cast<float>(n1), static_cast<float>(n2), theta, d};
                            diffs.push_back(diff);
                        };
                    };
                };
            };
        };
    };
    return diffs;
};

PYBIND11_MODULE(backend, m)
{
    m.doc() = "backend c++ implementation"; // optional module docstring
    py::bind_vector<std::vector<float>>(m, "FloatVector");
    py::bind_vector<std::vector<std::vector<float>>>(m, "FloatVector2D");
    m.def("find_coincidence", &find_coincidence, "A function that finds coincidence points.");
}
