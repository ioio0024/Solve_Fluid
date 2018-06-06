#ifndef MORTON_ORDER_H
#define MORTON_ORDER_H

#include "global_definitions.h"

#include <tuple>

inline Int64 z_encoding_impl(Int32 m) {
  Int64 n = m & 0x1fffff;
  n = (n | n << 32) & 0x1f00000000ffff;  // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
  n = (n | n << 16) & 0x1f0000ff0000ff;  // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
  n = (n | n <<  8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
  n = (n | n <<  4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
  n = (n | n <<  2) & 0x1249249249249249;
  return n;
}

inline Int64 z_encoding(Int32 i, Int32 j, Int32 k) {
  Int64 n = z_encoding_impl(i)
             | (z_encoding_impl(j) << 1)
             | (z_encoding_impl(k) << 2);
  return n;
}

inline Int32 z_decoding_impl(Int64 n) {
  n &= 0x1249249249249249;
  n = (n | n >> 2) & 0x10c30c30c30c30c3;
  n = (n | n >> 4) & 0x100f00f00f00f00f;
  n = (n | n >> 8) & 0x1f0000ff0000ff;
  n = (n | n >> 16) & 0x1f00000000ffff;
  n = (n | n >> 32) & 0x1fffff;
  return static_cast<Int32>(n);
}

inline std::tuple<Int32, Int32, Int32> z_decoding(Int64 n) {
  return std::make_tuple(z_decoding_impl(n),
                         z_decoding_impl(n >> 1),
                         z_decoding_impl(n >> 2));
}

#endif // MORTON_ORDER_H
