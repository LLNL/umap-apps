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

using pixel_type = float;
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

std::pair<double, std::vector<std::pair<pixel_type, vector_t>>>
shoot_vector(const median::cube_t<pixel_type> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::pair<pixel_type, vector_t>> result(num_random_vector);

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
    std::discrete_distribution<int> plus_or_minus{-1, 1};

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
      double x_slope = x_beta_dist(rnd_engine) * plus_or_minus(rnd_engine) * 2;
      double y_slope = y_beta_dist(rnd_engine) * plus_or_minus(rnd_engine) * 2;

      vector_t vector{x_intercept, x_slope, y_intercept, y_slope};

      vector_iterator<pixel_type> begin(cube, vector, 0);
      vector_iterator<pixel_type> end(vector_iterator<pixel_type>::create_end(cube, vector));

      // median calculation using Torben algorithm
      const auto start = utility::elapsed_time_sec();
      result[i].first = torben(begin, end);
      total_execution_time += utility::elapsed_time_sec(start);
      result[i].second = vector;
    }
  }

  return std::make_pair(total_execution_time, result);
}

void print_top_median(const median::cube_t<pixel_type> &cube, std::vector<std::pair<pixel_type, vector_t>> &result) {
  // Sort the results by the descending order of median value
  std::partial_sort(result.begin(),
                    result.begin() + 10, // get only top 10 elements
                    result.end(),
                    [](const std::pair<pixel_type, vector_t> &lhd,
                       const std::pair<pixel_type, vector_t> &rhd) {
                      return (lhd.first > rhd.first);
                    });

  // Print out the top 10 median values and corresponding pixel values
  std::cout << "Top 10 median and corresponding pixel values (NaN values are not used in median calculation)"
            << std::endl;
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.precision(2);
  for (size_t i = 0; i < 10; ++i) {
    const pixel_type median = result[i].first;
    const vector_t vector = result[i].second;
    std::cout << "[" << i << "]"
              << "\n Median: " << median
              << "\n Vector: " << vector.x_slope << " " << vector.x_intercept
              << " " << vector.y_slope << " " << vector.y_intercept << std::endl;

    vector_iterator<pixel_type> iterator(cube, vector, 0);
    vector_iterator<pixel_type> end(vector_iterator<pixel_type>::create_end(cube, vector));
    for (; iterator != end; ++iterator) {
        std::cout << " " << *iterator;
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

  size_t BytesPerElement;
  median::cube_t<pixel_type> cube;

  // Alloc memory space and read data from fits files with umap
  cube.data = (pixel_type *)utility::umap_fits_file::PerFits_alloc_cube(
      options.filename, &BytesPerElement, &cube.size_x, &cube.size_y, &cube.size_k);

  if (cube.data == nullptr) {
    std::cerr << "Failed to allocate memory for cube\n";
    std::abort();
  }

  if (sizeof(pixel_type) != BytesPerElement) {
    std::cerr << "Wrong pixel type" << std::endl;
    std::abort();
  }

  {
    const char *time_stamp_file_name = std::getenv("TIME_STAMP_FILE");
    if (time_stamp_file_name != nullptr) {
      std::ifstream ifs(time_stamp_file_name);
      if (!ifs.is_open()) {
        std::cerr << "Cannot open " << time_stamp_file_name << std::endl;
        std::abort();
      }
      std::vector<double> time_stamps;
      for (double time_stamp; ifs >> time_stamp;) {
        time_stamps.emplace_back(time_stamp);
      }
      if (time_stamps.size() != cube.size_k) {
        std::cerr << "#of lines in " << time_stamp_file_name << " is not the same as #of fits files" << std::endl;
        std::abort();
      }
      cube.time_stamps = new double[time_stamps.size()];
      std::memcpy(cube.time_stamps, time_stamps.data(), sizeof(double) * time_stamps.size());
    } else {
      // If a list of time stamps is not given, assume that the difference between two frames is 1.0
      cube.time_stamps = new double[cube.size_k];
      for (size_t i = 0; i < cube.size_k; ++i) cube.time_stamps[i] = i * 1.0;
    }
  }

  std::size_t num_random_vector = default_num_random_vector;
  {
    const char *buf = std::getenv("NUM_VECTORS");
    if (buf != nullptr) {
      num_random_vector = std::stoll(buf);
    }
  }

  auto result = shoot_vector(cube, num_random_vector);
  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << result.first
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / result.first << std::endl;
  print_top_median(cube, result.second);

  delete[] cube.time_stamps;
  utility::umap_fits_file::PerFits_free_cube(cube.data);

  return 0;
}
