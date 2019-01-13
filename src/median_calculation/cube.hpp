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


#ifndef UMAP_APPS_MEDIAN_CALCULATION_CUBE_HPP
#define UMAP_APPS_MEDIAN_CALCULATION_CUBE_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <tuple>

#include "utility.hpp"

#define MEDIAN_CALCULATION_COLUMN_MAJOR 1
#define MEDIAN_CALCULATION_VERBOSE_OUT_OF_RANGE 0

namespace median {

template <typename _pixel_type>
class cube {
 public:

  using pixel_type = _pixel_type;

  /// -------------------------------------------------------------------------------- ///
  /// Constructor
  /// -------------------------------------------------------------------------------- ///
  cube() = default;

  cube(const size_t size_x,
       const size_t size_y,
       const size_t size_k,
       pixel_type *const image_data,
       std::vector<double> timestamp_list)
      : m_size_x(size_x),
        m_size_y(size_y),
        m_size_k(size_k),
        m_image_data(image_data),
        m_timestamp_list(std::move(timestamp_list)) {
    assert(m_size_k <= m_timestamp_list.size());
  }

  ~cube() = default; // Default destructor
  cube(const cube &) = default; // Default copy constructor
  cube(cube &&) noexcept = default; // Default move constructor
  cube &operator=(const cube &) = default; // Default copy assignment
  cube &operator=(cube &&) = default; // Default move assignment

  /// -------------------------------------------------------------------------------- ///
  /// Public methods
  /// -------------------------------------------------------------------------------- ///

  /// \brief Return frame size
  size_t frame_size() const {
    return m_size_x * m_size_y;
  }

  /// \brief Return cube size
  size_t cube_size() const {
    return m_size_x * m_size_y * m_size_k;
  }

  /// \brief Returns TRUE if the given x-y-k coordinate points at out-of-range
  bool out_of_range(const ssize_t x, const ssize_t y, const ssize_t k) const {
    return (index_in_cube(x, y, k) == -1);
  }

  /// \brief Returns a pixel value of the given x-y-k coordinate
  /// A returned value can be NaN value
  pixel_type get_pixel_value(const ssize_t x, const ssize_t y, const ssize_t k) const {
    assert(!out_of_range(x, y, k));
    return reverse_byte_order<pixel_type>(m_image_data[index_in_cube(x, y, k)]);
  }

  /// \brief Returns the size of cube (x, y, k) in tuple
  std::tuple<size_t, size_t, size_t> size() const {
    return std::make_tuple(m_size_x, m_size_y, m_size_k);
  }

  const pixel_type* image_data() const {
    return m_image_data;
  }

  double timestamp(const size_t k) const {
    assert(k < m_timestamp_list.size());
    return m_timestamp_list[k];
  }

 private:
  /// -------------------------------------------------------------------------------- ///
  /// Private methods
  /// -------------------------------------------------------------------------------- ///

  /// \brief Returns the index of the given (x, y, k) in a 3D cube
  /// Returns -1 if the given xyz coordinate points at out side of the cube
  ssize_t index_in_cube(const ssize_t x, const ssize_t y, const ssize_t k) const {
    if (x < 0 || m_size_x <= x || y < 0 || m_size_y <= y) {
#if MEDIAN_CALCULATION_VERBOSE_OUT_OF_RANGE
      std::cerr << "Frame index is out-of-range: "
              << "(" << x << ", " << y << ") is out of "
              << "(" << m_size_x << ", " << m_size_y << ")" << std::endl;
#endif
      return -1;
    }

    if (k < 0 || m_size_k <= k) {
#if MEDIAN_CALCULATION_VERBOSE_OUT_OF_RANGE
      std::cerr << "Cube index is out-of-range: "
              << "(" << x << ", " << y << ", " << k << ") is out of "
              << "(" << m_size_x << ", " << m_size_y << ", " << m_size_k << ")" << std::endl;
#endif
      return -1;
    }

#if MEDIAN_CALCULATION_COLUMN_MAJOR
    const ssize_t frame_index = x + y * m_size_x;
#else
    const ssize_t frame_index = x * m_size_y + y;
#endif

    const ssize_t cube_index = frame_index + k * frame_size();

    return cube_index;
  }


  /// -------------------------------------------------------------------------------- ///
  /// Private fields
  /// -------------------------------------------------------------------------------- ///
  size_t m_size_x;
  size_t m_size_y;
  size_t m_size_k;

  pixel_type *const m_image_data;

  std::vector<double> m_timestamp_list; // an array of the timestamp of each frame.
};

} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_CUBE_HPP
