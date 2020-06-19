#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <iostream>
#include "logging.h"

using float_2d_vec = std::vector<std::vector<float>>;

float coincidence(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, int &m1, int &m2, int &n1, int &n2, float &theta)
{
    float v11 = a1[0] * m1 + a2[0] * m2;
    float v12 = a1[1] * m1 + a2[1] * m2;
    float v21 = std::cos(theta) * (b1[0] * n1 + b2[0] * n2) - std::sin(theta) * (b1[1] * n1 + b2[1] * n2);
    float v22 = std::sin(theta) * (b1[0] * n1 + b2[0] * n2) + std::cos(theta) * (b1[1] * n1 + b2[1] * n2);
    float norm = std::sqrt((v11 - v21) * (v11 - v21) + (v12 - v22) * (v12 - v22));
    return norm;
};

float_2d_vec find_coincidence(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, std::vector<float> angles, int Ntrans, float crit)
{
    float_2d_vec diffs;
    for (int i = -Ntrans; i < (Ntrans + 1); i++)
    {
        for (int j = -Ntrans; j < (Ntrans + 1); j++)
        {
            for (int k = -Ntrans; k < (Ntrans + 1); k++)
            {
                for (int l = -Ntrans; l < (Ntrans + 1); l++)
                {
                    for (float theta : angles)
                    {
                        float d = coincidence(a1, a2, b1, b2, i, j, k, l, theta);
                        bool iszero = ((i == 0) && (j == 0)) || ((k == 0) && l == 0);
                        if ((d < crit) && !iszero)
                        {
                            std::vector<float> diff{static_cast<float>(i), static_cast<float>(j), static_cast<float>(k), static_cast<float>(l), theta, d};
                            diffs.push_back(diff);
                        };
                    };
                };
            };
        };
    };
    return diffs;
};

std::vector<float> get_unique_angles(float_2d_vec &vec)
{
    std::vector<float> angles;
    for (float_2d_vec::size_type i = 0; i < vec.size(); ++i)
    {
        angles.push_back(vec[i][4]);
    };
    std::unordered_set<float> set;
    std::copy(angles.begin(), angles.end(), std::inserter(set, set.end()));
    angles.clear();
    for (auto &i : set)
    {
        angles.push_back(i);
    };
    std::sort(angles.begin(), angles.end());
    return angles;
};

float_2d_vec filter_2dvec_by_angle(float_2d_vec &results, float &angle)
{
    float_2d_vec farr;
    for (float_2d_vec::size_type i = 0; i < results.size(); i++)
    {
        float diff = std::abs(results[i][4] - angle);
        if (diff < 1e-10)
        {
            farr.push_back(results[i]);
        };
    };
    return farr;
};

float_2d_vec build_unique_pairs(float_2d_vec &farr, float &theta)
{
    float_2d_vec unpairs;
    for (float_2d_vec::size_type i = 0; i < farr.size(); ++i)
    {
        for (float_2d_vec::size_type j = 0; j < farr.size(); ++j)
        {
            if (j != i)
            {
                float m1 = farr[i][0];
                float m2 = farr[i][1];
                float n1 = farr[i][2];
                float n2 = farr[i][3];
                float m3 = farr[j][0];
                float m4 = farr[j][1];
                float n3 = farr[j][2];
                float n4 = farr[j][3];
                float detM, detN;
                detM = m1 * m4 - m2 * m3;
                detN = n1 * n4 - n2 * n3;
                if ((detM > 0) && (detN > 0))
                {
                    std::vector<float> subvec{theta, m1, m2, m3, m4, n1, n2, n3, n4};
                    unpairs.push_back(subvec);
                };
            };
        };
    };
    if (unpairs.size() > 0)
    {
        return unpairs;
    }
    else
    {
        return {};
    }
}

float_2d_vec build_all_unique_pairs(float_2d_vec &results, std::vector<float> &unangles)
{
    float_2d_vec pairs;
    float_2d_vec subpairs;
    float_2d_vec farr;
    for (auto theta : unangles)
    {
        farr = filter_2dvec_by_angle(results, theta);
        subpairs = build_unique_pairs(farr, theta);
        for (float_2d_vec::size_type j = 0; j < subpairs.size(); j++)
        {
            pairs.push_back(subpairs[j]);
        }
        subpairs.clear();
    }
    if (!pairs.empty())
    {
        return pairs;
    }
    else
    {
        return {};
    }
};

float angle_between_vectors(float theta, float &m1, float &m2, float &m3, float &m4, std::vector<float> &a1, std::vector<float> &a2)
{
    float c = std::cos(theta);
    float s = std::cos(theta);
    float v11 = m1 * (c * a1[0] - s * a1[1]) + m2 * (c * a2[0] - s * a2[1]);
    float v12 = m3 * (c * a1[0] - s * a1[1]) + m4 * (c * a2[0] - s * a2[1]);
    float w11 = m1 * (c * a1[1] + s * a1[0]) + m2 * (c * a2[1] + s * a2[0]);
    float w12 = m3 * (c * a1[1] + s * a1[0]) + m4 * (c * a2[1] + s * a2[0]);
    float normv = std::sqrt(v11 * v11 + v12 * v12);
    float normw = std::sqrt(w11 * w11 + w12 * w12);
    float dot = (v11 / normv) * (w11 / normw) + (v12 / normv) * (w12 / normw);
    float gamma = std::acos(dot) * 180 / 3.141592653589793238462643383279502884;
    return gamma;
};

// Function to return gcd of a and b
int gcd(int a, int b)
{
    a = std::abs(a);
    b = std::abs(b);
    int temp;
    while (b > 0)
    {
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
};

// Function to find gcd of a sorted vector
int find_gcd(std::vector<float> &subvec)
{
    int result = static_cast<int>(subvec[0]);
    for (std::vector<float>::size_type i = 0; i < subvec.size(); i++)
    {
        result = gcd(static_cast<int>(subvec[i]), result);
        if (result == 1)
        {
            return 1;
        }
    }
    return result;
};

float_2d_vec filter_unique_pairs_by_angle(float_2d_vec &pairs, std::vector<float> &a1, std::vector<float> &a2, std::vector<float> &b1, std::vector<float> &b2, float mingamma, float maxgamma)
{
    float_2d_vec fpairs;
    for (float_2d_vec::size_type i = 0; i < pairs.size(); i++)
    {
        // theta, m1, m2, m3, m4, n1, n2, n3, n4
        float theta = pairs[i][0];
        float m1 = pairs[i][1];
        float m2 = pairs[i][2];
        float m3 = pairs[i][3];
        float m4 = pairs[i][4];
        float n1 = pairs[i][5];
        float n2 = pairs[i][6];
        float n3 = pairs[i][7];
        float n4 = pairs[i][8];
        float gamma1 = angle_between_vectors(0, m1, m2, m3, m4, a1, a2); // this one is not rotated
        float gamma2 = angle_between_vectors(theta, n1, n2, n3, n4, b1, b2);
        if (((gamma1 >= mingamma) & (gamma1 <= maxgamma)) & ((gamma2 >= mingamma) & (gamma2 <= maxgamma)))
        {
            fpairs.push_back(pairs[i]);
        }
    }
    return fpairs;
};

float_2d_vec filter_unique_pairs_by_gcd(float_2d_vec &pairs)
{
    float_2d_vec fpairs;
    std::vector<float> subvec;
    int v_gcd;
    for (float_2d_vec::size_type i = 0; i < pairs.size(); i++)
    {
        // theta, m1, m2, m3, m4, n1, n2, n3, n4
        subvec = pairs[i];
        subvec.erase(subvec.begin()); // remove theta
        std::sort(subvec.begin(), subvec.end());
        v_gcd = find_gcd(subvec);
        if (v_gcd == 1)
        {
            fpairs.push_back(pairs[i]);
        }
    }
    return fpairs;
};

float_2d_vec backend_routine(std::vector<float> a1, std::vector<float> a2, std::vector<float> b1, std::vector<float> b2, std::vector<float> angles, int Ntrans, float crit, float mingamma, float maxgamma)
{
    float_2d_vec results;
    float_2d_vec pairs;
    std::vector<float> unangles;
    results = find_coincidence(a1, a2, b1, b2, angles, Ntrans, crit);
    std::cout << "INFO    │ Found " << results.size() << " matching lattice points." << std::endl;
    if (results.size() > 1)
    {
        unangles = get_unique_angles(results);
        pairs = build_all_unique_pairs(results, unangles);
        pairs = filter_unique_pairs_by_angle(pairs, a1, a2, b1, b2, mingamma, maxgamma);
        pairs = filter_unique_pairs_by_gcd(pairs);
        if (!pairs.empty())
        {
            std::cout << "INFO    │ Constructed " << pairs.size() << " linearly independent pairs." << std::endl;
        }
        else
        {
            std::cout << "ERROR   │ Found no linearly independent solutions." << std::endl;
        }
    };
    return pairs;
}
