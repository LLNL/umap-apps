/*
This file is part of UMAP.  For copyright information see the COPYRIGHT
file in the top level directory, or at
https://github.com/LLNL/umap/blob/master/COPYRIGHT
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License (as published by the Free
Software Foundation) version 2.1 dated February 1999.  This program is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the terms and conditions of the GNU Lesser General Public License
for more details.  You should have received a copy of the GNU Lesser General
Public License along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

/// Memo about median calculation
/// x: horizontal dimension of a frame at a certain time point
/// y: vertical dimension of a frame at a certain time point
/// k: time dimension
/// A cube is a set of 'k' frames

#ifndef UMAP_APPS_MEDIAN_CALCULATION_UTILITY_HPP
#define UMAP_APPS_MEDIAN_CALCULATION_UTILITY_HPP

#include <iostream>
#include <cmath>
#include <cassert>

namespace median {

/// \brief Reverses byte order
/// \tparam T Type of value; currently only 4 Byte types are supported
/// \param x Input value
/// \return Given value being reversed byte order
template <typename T>
T reverse_byte_order(const T x) {
  static_assert(sizeof(T) == 4, "T is not a 4 byte type");
  T reversed_x;
  const char *const p1 = reinterpret_cast<const char *>(&x);
  char *const p2 = reinterpret_cast<char *>(&reversed_x);
  p2[0] = p1[3];
  p2[1] = p1[2];
  p2[2] = p1[1];
  p2[3] = p1[0];

  return reversed_x;
}

template <typename pixel_type>
bool is_nan(const pixel_type value) {
  return std::isnan(value);
}

} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_UTILITY_HPP
