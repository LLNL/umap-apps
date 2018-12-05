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
#include "torben.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include "../utility/time.hpp"

using pixel_t = float;
constexpr size_t default_num_random_vector = 100000;

class beta_distribution {
 public:
  beta_distribution(double a, double b)
      : m_x_gamma(a, 1.0),
        m_y_gamma(b, 1.0) {}

  template <typename rnd_engine>
  double operator()(rnd_engine &engine) {
    double x = m_x_gamma(engine);
    double y = m_y_gamma(engine);
    return x / (x + y);
  }

 private:
  std::gamma_distribution<> m_x_gamma;
  std::gamma_distribution<> m_y_gamma;
};

void map_fits(const std::string& filename, median::cube_t<pixel_t> &cube) {
  size_t byte_per_element;
  // Map FITS files using UMap
  cube.data = (pixel_t *)utility::umap_fits_file::PerFits_alloc_cube(filename, &byte_per_element, &cube.size_x, &cube.size_y, &cube.size_k);

  if (cube.data == nullptr) {
    std::cerr << "Failed to allocate memory for cube" << std::endl;
    std::abort();
  }

  if (sizeof(pixel_t) != byte_per_element) {
    std::cerr << "Wrong pixel type" << std::endl;
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

void set_timestamp(median::cube_t<pixel_t> &cube) {
  const char *timestamp_file_name = std::getenv("TIMESTAMP_FILE");
  if (timestamp_file_name != nullptr) {
    std::ifstream ifs(timestamp_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open " << timestamp_file_name << std::endl;
      std::abort();
    }
    std::vector<double> timestamps;
    for (double timestamp; ifs >> timestamp;) {
      timestamps.emplace_back(timestamp);
    }
    if (timestamps.size() != cube.size_k) {
      std::cerr << "#of lines in " << timestamp_file_name << " is not the same as #of fits files" << std::endl;
      std::abort();
    }
    cube.timestamps = new double[timestamps.size()];
    std::memcpy(cube.timestamps, timestamps.data(), sizeof(double) * timestamps.size());
  } else {
    // If a list of timestamps is not given, assume that the difference between two frames is 1.0
    cube.timestamps = new double[cube.size_k];
    for (size_t i = 0; i < cube.size_k; ++i) cube.timestamps[i] = i * 1.0;
  }
}

std::pair<double, std::vector<std::pair<pixel_t, vector_t>>>
shoot_vector(const median::cube_t<pixel_t> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::pair<pixel_t, vector_t>> result(num_random_vector);

  double total_execution_time = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
    std::mt19937 rnd_engine(123 + omp_get_thread_num());
#else
    std::mt19937 rnd_engine(123);
#endif
    std::uniform_int_distribution<int> x_start_dist(0, cube.size_x - 1);
    std::uniform_int_distribution<int> y_start_dist(0, cube.size_y - 1);
    beta_distribution x_beta_dist(3, 2);
    beta_distribution y_beta_dist(3, 2);
    std::uniform_int_distribution<int> plus_or_minus(0, 1);

    // Shoot random vectors using multiple threads
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < num_random_vector; ++i) {
      double x_intercept = x_start_dist(rnd_engine);
      double y_intercept = y_start_dist(rnd_engine);

      // Changed to the const value to 2 from 25 so that vectors won't access
      // out of range of the cube with a large number of frames
      //
      // This is a temporary measures
      double x_slope = x_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);
      double y_slope = y_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);

      vector_t vector{x_intercept, x_slope, y_intercept, y_slope};

      vector_iterator<pixel_t> begin(cube, vector, 0);
      vector_iterator<pixel_t> end(vector_iterator<pixel_t>::create_end(cube, vector));

      // median calculation using Torben algorithm
      const auto start = utility::elapsed_time_sec();
      result[i].first = torben(begin, end);
      total_execution_time += utility::elapsed_time_sec(start);
      result[i].second = vector;
    }
  }

  return std::make_pair(total_execution_time, result);
}

void print_top_median(const median::cube_t<pixel_t> &cube,
                      const size_t num_top,
                      std::vector<std::pair<pixel_t, vector_t>> &result) {
  // Sort the results by the descending order of median value
  std::partial_sort(result.begin(),
                    result.begin() + num_top, // get only top 10 elements
                    result.end(),
                    [](const std::pair<pixel_t, vector_t> &lhd,
                       const std::pair<pixel_t, vector_t> &rhd) {
                      return (lhd.first > rhd.first);
                    });

  // Print out the top 'num_top' median values and corresponding pixel values
  std::cout << "Top " << num_top << " median and pixel values (skip NaN value)" << std::endl;
  for (size_t i = 0; i < num_top; ++i) {
    const pixel_t median = result[i].first;
    const vector_t vector = result[i].second;
    std::cout << "[" << i << "]" << std::endl;
    std::cout << "Median: " << median << std::endl;
    std::cout << "Vector (x-slope, x-intercept, y-slope, y-intercept): "
              << vector.x_slope << ", " << vector.x_intercept << ", " << vector.y_slope << ", " << vector.y_intercept << std::endl;

    std::cout << "Values (x, y, k):" << std::endl;
    vector_iterator<pixel_t> iterator(cube, vector, 0);
    vector_iterator<pixel_t> end(vector_iterator<pixel_t>::create_end(cube, vector));
    for (; iterator != end; ++iterator) {
      ssize_t x, y, k;
      std::tie(x, y, k) = iterator.coordinate();
      std::cout << " [ " << x << ", " << y << ", " << k << " ] = " << *iterator << std::endl;
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

  median::cube_t<pixel_t> cube;
  map_fits(options.filename, cube);
  set_timestamp(cube);
  const std::size_t num_random_vector = get_num_vectors();

  auto result = shoot_vector(cube, num_random_vector);
  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << result.first
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / result.first << std::endl;

  print_top_median(cube, std::min(num_random_vector, static_cast<size_t>(10)), result.second);

  delete[] cube.timestamps;
  utility::umap_fits_file::PerFits_free_cube(cube.data);

  return 0;
}
