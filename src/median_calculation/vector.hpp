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

#ifndef UMAP_APPS_MEDIAN_CALCULATION_VECTOR_HPP
#define UMAP_APPS_MEDIAN_CALCULATION_VECTOR_HPP

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>
#include <unordered_set>

#include <umap/umap.h>

#include "utility.hpp"
#include "cube.hpp"

namespace median {

template <typename pixel_type>
struct vector_xy {
  double x_slope;
  uint64_t x_intercept;
  double y_slope;
  uint64_t y_intercept;
  double z_intercept;
  pixel_type* pixels;
  uint64_t npixels;
  vector_xy() = default;
  vector_xy(const vector_xy&) = default;
  vector_xy(double x_slope, uint64_t x_intercept, double y_slope, uint64_t y_intercept, uint64_t z_intercept, pixel_type* pixels)
    : x_slope(x_slope), y_slope(y_slope), x_intercept(x_intercept), y_intercept(y_intercept), z_intercept(z_intercept), pixels(pixels) {}

  void set(double x_slope, uint64_t x_intercept, double y_slope, uint64_t y_intercept, uint64_t z_intercept, pixel_type* pixels) {
    this->x_slope = x_slope;
    this->x_intercept = x_intercept;
    this->y_slope = y_slope;
    this->y_intercept = y_intercept;
    this->z_intercept = z_intercept;
    this->pixels = pixels;
  }
  /// \brief Returns the xy position at a given offset
  std::pair<uint64_t, uint64_t> position(const double offset) const {
    const uint64_t x = std::round(x_slope * (z_intercept-offset) + x_intercept);
    const uint64_t y = std::round(y_slope * (z_intercept-offset) + y_intercept);
    return std::make_pair(x, y);
  }
};

} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_VECTOR_HPP

// -------------------------------------------------------------------------------- //
// Example of a function that returns pixel value by a window.
// To use this function, just call it in '*' operator
// value_type operator*() const {
//   return get_pixel_value_with_window()
// }
// -------------------------------------------------------------------------------- //
/*
pixel_type get_pixel_value_with_window(const uint64_t window_size_x, const uint64_t window_size_y) const {
  const auto xy = current_xy_position();
  uint64_t x = xy.first - window_size_x / 2;
  uint64_t y = xy.second - window_size_y / 2;

  pixel_type result = 0;
  uint64_t num_valid_values = 0;
  for (uint64_t offset_y = 0; offset_y < window_size_y; ++offset_y) {
    for (uint64_t offset_x = 0; offset_x < window_size_x; ++offset_x) {
      if (m_cube.out_of_range(x + offset_x, y + offset_y, m_current_k_pos)) continue;

      const pixel_type value = m_cube.get_pixel_value(x + offset_x, y + offset_y, m_current_k_pos);
      if (is_nan(value)) continue;

      result += value;
      ++num_valid_values;
    }
  }

  return (double)result / num_valid_values;
}
*/
