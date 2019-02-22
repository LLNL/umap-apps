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
#include <math.h>

#include "utility.hpp"
#include "cube.hpp"

namespace median {

struct vector_xy {
  double x_slope;
  double x_intercept;
  double y_slope;
  double y_intercept;

  /// \brief Returns the xy position at a given offset
  std::pair<ssize_t, ssize_t> position(const double offset) const {
    const ssize_t x = std::round(x_slope * offset + x_intercept);
    const ssize_t y = std::round(y_slope * offset + y_intercept);
    return std::make_pair(x, y);
  }
};

// Iterator class to use the Torben function with vector model
// This class is a minimum implementation of an iterator to use the Torben function
template <typename pixel_type>
class cube_iterator_with_vector {
 public:
  using value_type = pixel_type;

  /// -------------------------------------------------------------------------------- ///
  /// Constructors
  /// -------------------------------------------------------------------------------- ///

  // Configured as an iterator pointing to the 'end'
  cube_iterator_with_vector(const cube<pixel_type> &_cube,
                            const vector_xy &_vector_xy)
      : m_cube(_cube),
        m_vector(_vector_xy),
        m_current_k_pos(std::get<2>(m_cube.size())) {}

  cube_iterator_with_vector(cube<pixel_type> _cube,
                            vector_xy _vector_xy,
                            size_t _start_k_pos)
      : m_cube(std::move(_cube)),
        m_vector(_vector_xy),
        m_current_k_pos(_start_k_pos) {

    // m_current_k_pos must be less than size_k always
    const size_t size_k = std::get<2>(m_cube.size());
    if (size_k < m_current_k_pos) {
      m_current_k_pos = size_k;
      return;
    }

    // Move to the first valid pixel
    const auto xy = current_xy_position();
    if (m_cube.out_of_range(xy.first, xy.second, m_current_k_pos) // This one has to be evaluated first
        || is_nan(m_cube.get_pixel_value(xy.first, xy.second, m_current_k_pos))) {
      move_to_next_valid_pixel();
    }
  }

  // Use default copy constructor
  cube_iterator_with_vector(const cube_iterator_with_vector &) = default;

  /// -------------------------------------------------------------------------------- ///
  /// Operators and public methods
  /// -------------------------------------------------------------------------------- ///

  // To support
  // iterator1 == iterator2
  bool operator==(const cube_iterator_with_vector &other) const {
    return m_current_k_pos == other.m_current_k_pos;
  }

  // To support
  // iterator1 != iterator2
  bool operator!=(const cube_iterator_with_vector &other) const {
    return !(*this == other);
  }

  // To support
  // value_type val = *iterator
  // Note: it is possible for this to return nan and functions dependent on this should take necessary precautions
  value_type operator*() const {
    return get_pixel_value_with_streak();
  }

  // To support
  // ++iterator
  cube_iterator_with_vector &operator++() {
    move_to_next_valid_pixel();
    return (*this);
  }

 private:
  /// -------------------------------------------------------------------------------- ///
  /// Private methods
  /// -------------------------------------------------------------------------------- ///
  std::pair<ssize_t, ssize_t> current_xy_position() const {
    const double time_offset = m_cube.timestamp(m_current_k_pos) - m_cube.timestamp(0);
    return m_vector.position(time_offset);
  }
  
  pixel_type get_pixel_value_with_streak() const {
    const auto xy = current_xy_position();
    ssize_t x_pos = xy.first;
    ssize_t y_pos = xy.second;

    double exp_time = m_cube.exposuretime(m_current_k_pos);
    size_t psf_width = m_cube.psf(m_current_k_pos);

    size_t streak_length = sqrt(pow((exp_time / m_vector.x_slope), 2) + pow((exp_time / m_vector.y_slope), 2));
    double phi = atan2(m_vector.x_slope,m_vector.y_slope); ///simplified from (streak_length_y/streak_length_x)

    pixel_type result = 0;

    for (size_t x_offset = floor(-streak_length/2 -psf_width); x_offset <= ceil(streak_length/2 + psf_width); ++x_offset) {
      for (size_t y_offset = floor(-psf_width); y_offset <= ceil(psf_width); ++y_offset) {
        if (m_cube.out_of_range(x_pos + x_offset, y_pos + y_offset, m_current_k_pos)) continue;

        size_t x_pixel = std::round(x_pos + cos(phi)*x_offset - sin(phi)*y_offset);
        size_t y_pixel = std::round(y_pos + sin(phi)*x_offset + cos(phi)*y_offset);
        
		const pixel_type value = m_cube.get_pixel_value(x_pixel, y_pixel, m_current_k_pos);
        if (is_nan(value)) continue;
        
		///Following code is for a weighted sum using convolution of gaussian with streak
        pixel_type psfwidth_sq = (pixel_type) psf_width*psf_width;
        pixel_type coef = (1/(2*streak_length))*(2*psfwidth_sq*M_PI/sqrt(2*M_PI*psfwidth_sq));
        pixel_type xarg_denom = 2*sqrt(2*psfwidth_sq);
        pixel_type xarg1 = (streak_length - 2*x_offset)/xarg_denom;
        pixel_type xarg2 = (streak_length + 2*x_offset)/xarg_denom;
        pixel_type xterm = erf(xarg1) + erf(xarg2);
        pixel_type yterm = exp(-0.5 *y_offset*y_offset/psfwidth_sq);
        pixel_type weight = coef*xterm*yterm;

        result += value*weight;
    }
  }

  if (result == 0) result = nan; // if all streak pixels are nan, return nan
  return result;
}

  // Find the next non-NaN value
  void move_to_next_valid_pixel() {
    ++m_current_k_pos;
    const size_t size_k = std::get<2>(m_cube.size());

    for (; m_current_k_pos < size_k; ++m_current_k_pos) {
      const auto xy = current_xy_position();

      if (!m_cube.out_of_range(xy.first, xy.second, m_current_k_pos)) return;

	  //We're no longer skipping nan values here, but instead after the *() operator.
	  //The exact places where this is done are in the vector_sum and torben functions.
	  //This allows us to pull potentially relevant info from around the exact pixel 
	  //that a vector intersects if the exact pixel is nan.

      //if (!is_nan(m_cube.get_pixel_value(xy.first, xy.second, m_current_k_pos))) return;
    }

    m_current_k_pos = size_k; // Prevent the case, m_current_k_pos > size_k.
  }

  /// -------------------------------------------------------------------------------- ///
  /// Private fields
  /// -------------------------------------------------------------------------------- ///
  const cube<pixel_type> m_cube;
  vector_xy m_vector;
  size_t m_current_k_pos;
};

} // namespace median

#endif //UMAP_APPS_MEDIAN_CALCULATION_VECTOR_HPP
