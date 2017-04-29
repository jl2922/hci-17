#ifndef HCI_TYPES_H_
#define HCI_TYPES_H_

#include "std.h"

typedef std::int8_t TinyInt;  // [-128, 127].
typedef std::uint16_t Orbital;  // [0, 65535 (2^16-1)].
typedef std::uint64_t BigUnsignedInt;  // [0, >1.84e19 (2^64-1)].
typedef std::array<int, 3> Int3;
typedef std::pair<int, int> IntPair;
typedef std::array<TinyInt, 3> TinyInt3;
typedef std::pair<TinyInt3, double> TinyInt3Double;
typedef std::vector<int> Ints;
typedef std::vector<Orbital> Orbitals;

#endif