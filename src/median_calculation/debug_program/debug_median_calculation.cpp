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

/// \brief This program is designed to debug our median calculation code WITHOUT UMap

#include <iostream>
#include <string>

#include "../../utility/mmap.hpp"
#include "../torben.hpp"
#include "../utility.hpp"
#include "../vector.hpp"
#include "../cube.hpp"

using pixel_type = float;
using namespace median;

constexpr size_t size_x = 32;
constexpr size_t size_y = 32;
constexpr size_t size_k = 10;

const std::string file_name("debug_file");

pixel_type *map_data(const size_t size_x, const size_t size_y, const size_t size_k,
                     const std::string &file_name) {
  /// ----- create and map output (graph) file ----- ///
  const size_t file_size = size_x * size_y * size_k * sizeof(pixel_type);
  if (!utility::create_file(file_name)) {
    std::cerr << "Failed to create a file: " << file_name << std::endl;
    std::abort();
  }
  if (!utility::extend_file_size(file_name, file_size)) {
    std::cerr << "Failed to extend a file: " << file_name << std::endl;
    std::abort();
  }

  int fd = -1;
  void *map_raw_address = nullptr;
  std::tie(fd, map_raw_address) = utility::map_file_write_mode(file_name, nullptr, file_size, 0);
  if (fd == -1 || map_raw_address == nullptr) {
    std::cerr << "Failed to map a file with write mode: " << file_name << std::endl;
    std::abort();
  }
  ::close(fd);

  return static_cast<pixel_type *>(map_raw_address);
}

void init_data(const size_t size_x, const size_t size_y, const size_t size_k, pixel_type *data) {

  for (size_t x = 0; x < size_x; ++x) {
    for (size_t y = 0; y < size_y; ++y) {
      for (size_t k = 0; k < size_k; ++k) {
        const pixel_type value = static_cast<pixel_type>(x * 100 + y + k / 100.0);
        data[x + size_x * y + size_x * size_y * k] = reverse_byte_order(value);
      }
    }
  }
}

void calculate_median(const cube<pixel_type> &cube, const vector_xy &vector) {
  cube_iterator_with_vector<pixel_type> begin(cube, vector, 0);
  cube_iterator_with_vector<pixel_type> end(cube, vector);

  // median calculation w/ Torben algorithm
  const auto median_val = torben(begin, end);

  std::cout << "Median: " << median_val << std::endl;
  std::cout << "Vector (x-slope, x-intercept, y-slope, y-intercept): "
            << vector.x_slope << ", " << vector.x_intercept << ", " << vector.y_slope << ", " << vector.y_intercept
            << std::endl;

  std::cout << "[x, y, k] =\tPixel values" << std::endl;
  for (size_t k = 0; k < std::get<2>(cube.size()); ++k) {
    const double time_offset = cube.timestamp(k) - cube.timestamp(0);
    const ssize_t x = vector.position(time_offset).first;
    const ssize_t y = vector.position(time_offset).second;

    std::cout << " [ " << x << ", " << y << ", " << k << " ] =\t";

    if (cube.out_of_range(x, y, k)) std::cout << "OOR" << std::endl; // Out of Range
    else std::cout << cube.get_pixel_value(x, y, k) << std::endl;
  }
  std::cout << std::endl;
}

int main() {
  pixel_type *data = map_data(size_x, size_y, size_k, file_name);
  init_data(size_x, size_y, size_k, data);

  std::vector<double> timestamp_list(size_k);
  for (size_t i = 0; i < size_k; ++i) timestamp_list[i] = i * 1.0;

  cube<pixel_type> cube(size_x, size_y, size_k, data, std::move(timestamp_list));

  calculate_median(cube, vector_xy{1, 0, 1, 0});

  calculate_median(cube, vector_xy{2, 10, 10, 5});

  calculate_median(cube, vector_xy{2, -2, 10, 5});

  return 0;
}
