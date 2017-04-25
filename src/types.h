#ifndef HCI_TYPES_H_
#define HCI_TYPES_H_

#include <array>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>  // std::pair.
#include <vector>

using std::array;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::pair;
using std::vector;

typedef signed char TinyInt;  // [-128, 127].
typedef unsigned long long BigUnsignedInt;  // [0, 2^64-1 (~1.84e19)].
typedef array<int, 3> Int3;
typedef array<TinyInt, 3> TinyInt3;
typedef pair<TinyInt3, double> TinyInt3Double;

const double EPSILON = 1.0e-20;

#endif