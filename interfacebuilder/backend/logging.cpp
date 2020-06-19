#include <iostream>
#include <vector>
#include "logging.h"

void log_2dvec(std::vector<std::vector<float>> &vec)
{
    for (std::vector<std::vector<float>>::size_type i = 0; i < vec.size(); i++)
    {
        for (std::vector<float>::size_type j = 0; j < vec[i].size(); j++)
        {
            std::cout << vec[i][j] << ' ';
        }
        std::cout << std::endl;
    }
};

void log_1dvec(std::vector<float> &vec)
{
    for (std::vector<float>::size_type i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << std::endl;
    };
};