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
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utility/commandline.hpp"
#include "../utility/umap_fits_file.hpp"
#include "../utility/time.hpp"
#include "torben.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include "cube.hpp"
#include "beta_distribution.hpp"

using namespace median;

using pixel_type = float;
constexpr size_t default_num_random_vector = 100000;

void map_fits(const std::string &filename,
              size_t *size_x,
              size_t *size_y,
              size_t *size_k,
              pixel_type **image_data) {
  size_t byte_per_element;
  // Map FITS files using UMap
  *image_data = (pixel_type *)utility::umap_fits_file::PerFits_alloc_cube(filename, &byte_per_element,
                                                                          size_x, size_y, size_k);

  if (*image_data == nullptr) {
    std::cerr << "Failed to allocate memory for cube" << std::endl;
    std::abort();
  }

  if (sizeof(pixel_type) != byte_per_element) {
    std::cerr << "Pixel type is not float" << std::endl;
    std::abort();
  }
}

std::size_t get_num_vectors() {
  std::size_t num_random_vector = default_num_random_vector;
  const char *buf = std::getenv("NUM_VECTORS");
  if (buf != nullptr) {
    num_random_vector = std::stoll(buf);
  }
  return num_random_vector;
}

std::vector<double> read_timestamp(const size_t size_k) {
  std::vector<double> timestamp_list;

  const char *timestamp_file_name = std::getenv("TIMESTAMP_FILE");
  if (timestamp_file_name != nullptr) {
    std::ifstream ifs(timestamp_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open " << timestamp_file_name << std::endl;
      std::abort();
    }
    for (double timestamp; ifs >> timestamp;) {
      timestamp_list.emplace_back(timestamp);
    }
    if (timestamp_list.size() != size_k) {
      std::cerr << "#of lines in " << timestamp_file_name << " is not the same as #of fits files" << std::endl;
      std::abort();
    }
  } else {
    // If a list of timestamps is not given, assume that the difference between two frames is 1.0
    timestamp_list.resize(size_k);
    for (size_t i = 0; i < size_k; ++i) timestamp_list[i] = i * 1.0;
  }

  return timestamp_list;
}

std::pair<double, std::vector<std::pair<pixel_type, vector_xy>>>
shoot_vector(const cube<pixel_type> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::pair<pixel_type, vector_xy>> result(num_random_vector);

  double total_execution_time = 0.0;
  int numthreads = 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP

    if (omp_get_thread_num() == 0) {
      numthreads = omp_get_num_threads();
      std::cout << "num threads: " << numthreads << std::endl;
    }
    std::mt19937 rnd_engine(123 + omp_get_thread_num());
#else
    std::mt19937 rnd_engine(123);
#endif
    std::uniform_int_distribution<int> x_start_dist(0, std::get<0>(cube.size()) - 1);
    std::uniform_int_distribution<int> y_start_dist(0, std::get<1>(cube.size()) - 1);
    beta_distribution x_beta_dist(3, 2);
    beta_distribution y_beta_dist(3, 2);
    std::uniform_int_distribution<int> plus_or_minus(0, 1);

    // Shoot random vectors using multiple threads
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < num_random_vector; ++i) {
      const double x_intercept = x_start_dist(rnd_engine);
      const double y_intercept = y_start_dist(rnd_engine);

      // Changed to the const value to 2 from 25 so that vectors won't access
      // out of range of the cube with a large number of frames
      //
      // This is a temporary measures
      const double x_slope = x_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);
      const double y_slope = y_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);

      vector_xy vector{x_slope, x_intercept, y_slope, y_intercept};

      cube_iterator_with_vector<pixel_type> begin(cube, vector, 0.0);
      cube_iterator_with_vector<pixel_type> end(cube, vector);

      // median calculation using Torben algorithm
      const auto start = utility::elapsed_time_sec();
      result[i].first = torben(begin, end);
      total_execution_time += utility::elapsed_time_sec(start);
      result[i].second = vector;
    }
  }

  return std::make_pair(total_execution_time / numthreads, result);
}

void print_top_median(const cube<pixel_type> &cube,
                      const size_t num_top,
                      std::vector<std::pair<pixel_type, vector_xy>> &result) {

  // Sort the results by the descending order of median value
  std::sort(result.begin(), result.end(),
            [](const std::pair<pixel_type, vector_xy> &lhd,
               const std::pair<pixel_type, vector_xy> &rhd) {
              return (lhd.first > rhd.first);
            });

  // Print out the top 'num_top' median values and corresponding pixel values
  std::cout << "Top " << num_top << " median and pixel values (skip NaN value)" << std::endl;
  for (size_t i = 0; i < num_top; ++i) {
    const pixel_type median = result[i].first;
    const vector_xy vector = result[i].second;

    std::cout << "[" << i << "]" << std::endl;
    std::cout << "Median: " << median << std::endl;
    std::cout << "Vector (x-slope, x-intercept, y-slope, y-intercept): "
              << vector.x_slope << ", " << vector.x_intercept
              << ", " << vector.y_slope << ", " << vector.y_intercept << std::endl;

    std::cout << "Values (x, y, k):" << std::endl;
    for (size_t k = 0; k < std::get<2>(cube.size()); ++k) {
      const double time_offset = cube.timestamp(k) - cube.timestamp(0);
      const ssize_t x = vector.position(time_offset).first;
      const ssize_t y = vector.position(time_offset).second;

      std::cout << " [ " << x << ", " << y << ", " << k << " ] = ";

      if (cube.out_of_range(x, y, k)) std::cout << "OOR" << std::endl; // Out of Range
      else std::cout << cube.get_pixel_value(x, y, k) << std::endl;
    }
    std::cout << std::endl;
  }
}

int main(int argc, char **argv) {
  utility::umt_optstruct_t options;
  umt_getoptions(&options, argc, argv);

#ifdef _OPENMP
  omp_set_num_threads(options.numthreads);
#endif

  size_t size_x; size_t size_y; size_t size_k;
  pixel_type *image_data;
  map_fits(options.filename, &size_x, &size_y, &size_k, &image_data);
  cube<pixel_type> cube(size_x, size_y, size_k, image_data, read_timestamp(size_k));

  const std::size_t num_random_vector = get_num_vectors();

  const auto start = utility::elapsed_time_sec();
  auto result = shoot_vector(cube, num_random_vector);
  double txt = utility::elapsed_time_sec(start);
  double thread_exec = result.first;

  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << txt
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / txt << std::endl;

  print_top_median(cube, std::min(num_random_vector, static_cast<size_t>(10)), result.second);

  utility::umap_fits_file::PerFits_free_cube(image_data);

  return 0;
}
