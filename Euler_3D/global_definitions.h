#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <cstdint>
#include <vector>

using Int32 = std::int32_t;
using Int64 = std::int64_t;

using Real = double;
using Rvec = std::vector<Real>;

constexpr Int32 GHOST = 2;
constexpr Real CFL = 0.9;
constexpr Real EPS = 1.0e-8;

#endif // GLOBAL_DEFINITIONS_H
