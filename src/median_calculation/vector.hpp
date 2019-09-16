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

// Iterator class to use the Torben function with vector model
// This class is a minimum implementation of an iterator to use the Torben function
template <typename pixel_type>
class multi_vector {
 public:
  using value_type = pixel_type;

  /// -------------------------------------------------------------------------------- ///
  /// Constructors
  /// -------------------------------------------------------------------------------- ///

  multi_vector(cube<pixel_type> _cube, vector_xy<pixel_type>* _vector_xy,  uint64_t nvecs)
      : m_cube(std::move(_cube)),
        m_vector(std::move(_vector_xy)),
        num_vecs(nvecs),
        start_k(0) {
    page_size = umapcfg_get_umap_page_size();
    m_current_k_pos = start_k;
    vec_max = m_vector.size();
    populate_layers(start_k, std::get<2>(m_cube.size()));
  }

  // Use default copy constructor
  multi_vector(const multi_vector &) = default;

  /// -------------------------------------------------------------------------------- ///
  /// Operators and public methods
  /// -------------------------------------------------------------------------------- ///

  void fetch() {
    const uint64_t pf_size = 32;
    uint64_t i = 0, size_k = std::get<2>(m_cube.size());
    prefetch_layers(pf_size);
    for (i = pf_size; i < size_k - pf_size; i+=pf_size) {
      prefetch_layers(pf_size);
      populate_layers(i - pf_size, pf_size);
    }
    populate_layers(i, size_k - i);
  }

  uint64_t size() { return num_vecs; }

  const vector_xy<pixel_type>& get_vector(uint64_t i) const { return m_vector[i]; }
  const std::vector<vector_xy<pixel_type>>& get_vector() const { return m_vector; }

 private:
  /// -------------------------------------------------------------------------------- ///
  /// Private methods
  /// -------------------------------------------------------------------------------- ///

  void populate_layers(uint64_t start, uint64_t n) {
    const uint64_t size_k = std::get<2>(m_cube.size());
    for (uint64_t k = start; (k < start + n) && (k < size_k); k++) {
      const double time_offset = m_cube.timestamp(k) - m_cube.timestamp(0);
#pragma omp parallel for schedule(static, 64)
      for (uint64_t i = 0; i < vec_max; i++) {
        vector_xy<pixel_type>& v = m_vector[i];
        const auto xy = v.position(time_offset);
        if (time_offset >= v.z_intercept && !m_cube.out_of_range(xy.first, xy.second, k)) {
          v.pixels[v.npixels++] = *m_cube.get_pixel_addr(xy.first, xy.second, k);
        }
      }
      if (k % 100 == 0) { std::cout << k << " / " << size_k << std::endl; }
    }
  }

  void prefetch_layers(uint64_t n) {
    std::unordered_set<void*> addr_set;
    const uint64_t size_k = std::get<2>(m_cube.size());
    for (uint64_t k = m_current_k_pos; (k < m_current_k_pos + n) && (k < size_k); ++k) {
      const double time_offset = m_cube.timestamp(k) - m_cube.timestamp(0);
      for (uint64_t i = 0; i < vec_max; i++) {
        vector_xy<pixel_type>& v = m_vector[i];
        const auto xy = v.position(time_offset);
        if (m_cube.out_of_range(xy.first, xy.second, k)) {
          // Invalidate this vector
          // Move this vector to the end of the list and rerun the iteration
          vec_max -= 1; // mark one element at the end as invalid
          vector_xy<pixel_type> tmp = m_vector[i];
          m_vector[i] = m_vector[vec_max];
          m_vector[vec_max] = tmp;
          i -= 1;
          continue;
        }
        pixel_type* addr = m_cube.get_pixel_addr(xy.first, xy.second, k);
        addr_set.insert((void*)(((uintptr_t)addr) & ~(page_size - 1)));
      }
    }
    umap_prefetch_item* pf_list = (umap_prefetch_item*)malloc( sizeof(umap_prefetch_item) * addr_set.size() );
    assert(pf_list != 0);
    int i = 0;
    for (auto ptr : addr_set) {
      pf_list[i++].page_base_addr = ptr;
    }
    umap_prefetch(i, pf_list);
    free(pf_list);
    this->pf_k = m_current_k_pos + n;
  }

  /// -------------------------------------------------------------------------------- ///
  /// Private fields
  /// -------------------------------------------------------------------------------- ///
  const cube<pixel_type> m_cube;
  std::vector<vector_xy<pixel_type>> m_vector;
  uint64_t m_current_k_pos, start_k, pf_k, page_size, vec_max;
  uint64_t num_vecs;
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
