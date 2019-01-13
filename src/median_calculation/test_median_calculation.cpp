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

#include <iostream>
#include <cmath>

#include "../utility/commandline.hpp"
#include "../utility/umap_fits_file.hpp"
#include "torben.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include "cube.hpp"

using pixel_type = float;

constexpr size_t num_vectors = 6;
const pixel_type x_intercept[num_vectors] = {1058.2, 1325.606, 1010.564, 829.674, 1390.826, 1091.015};
const pixel_type x_slope[num_vectors] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
const pixel_type y_intercept[num_vectors] = {124, 424, 724, 1024, 1324, 1624};
const pixel_type y_slope[num_vectors] = {0, 0, 0, 0, 0, 0};
const pixel_type correct_median[num_vectors] = {14913.25, 15223.21, 2284.29, 8939.15, 24899.55, 2395.80};

using namespace median;

int main(int argc, char** argv)
{
  utility::umt_optstruct_t options;
  umt_getoptions(&options, argc, argv);

  size_t BytesPerElement;
  size_t size_x; size_t size_y; size_t size_k;
  pixel_type *image_data;

  image_data = (pixel_type*)utility::umap_fits_file::PerFits_alloc_cube(options.filename, &BytesPerElement, &size_x, &size_y, &size_k);

  std::vector<double> timestamp_list(size_k);
  for (size_t i = 0; i < size_k; ++i) timestamp_list[i] = i * 1.0;

  cube<pixel_type> cube(size_x, size_y, size_k, image_data, timestamp_list);

  for (int i = 0; i < num_vectors; ++i) {
    std::cout << "Vector " << i << std::endl;

    vector_xy vector{x_slope[i], x_intercept[i], y_slope[i], y_intercept[i]};
    cube_iterator_with_vector<pixel_type> begin(cube, vector, 0);
    cube_iterator_with_vector<pixel_type> end(cube, vector);

    // median calculation w/ Torben algorithm
    const auto median_val = torben(begin, end);

    // Check the result
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(2);
    if (std::fabs(median_val - correct_median[i]) < 0.01) {
      std::cout << " Correct " <<  median_val << " == " << correct_median[i] << std::endl;
    } else {
      std::cerr << " Error " <<  median_val << " != " << correct_median[i] << std::endl;

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

      std::abort();
    }
  }

  utility::umap_fits_file::PerFits_free_cube(image_data);

  std::cout << "Passed all tests" << std::endl;

  return 0;
}
