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
#include <iterator>

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

#define NGEN 4096
#define NPF 256

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

std::pair<double, std::vector<std::pair<pixel_type, vector_xy<pixel_type>>>>
shoot_vector(const cube<pixel_type> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::pair<pixel_type, vector_xy<pixel_type>>> result(num_random_vector);

  double total_execution_time = 0.0;
  int numthreads = 1;

#pragma omp parallel
  {
    if (omp_get_thread_num() == 0) {
      numthreads = omp_get_num_threads();
      std::cout << "num threads: " << numthreads << std::endl;
    }
  }
  const uint64_t batch_sz = num_random_vector;

  uint64_t rnd = 0;
  asm volatile("rdrand %0" : "=r"(rnd) ::);
  std::mt19937 rnd_engine(rnd);
  uint64_t size_x = std::get<0>(cube.size()), size_y = std::get<1>(cube.size()), size_k = std::get<2>(cube.size());
  std::uniform_int_distribution<uint64_t> x_start_dist(0, size_x - 1);
  std::uniform_int_distribution<uint64_t> y_start_dist(0, size_y - 1);
  std::uniform_real_distribution<double> z_start_dist(cube.timestamp(0), cube.timestamp(size_k - 1));
  beta_distribution x_beta_dist(3, 2);
  beta_distribution y_beta_dist(3, 2);
  std::uniform_int_distribution<int> plus_or_minus(0, 1);

  uint64_t pg_per_vec = size_k * sizeof(pixel_type) / 4096;
  uint64_t sz_per_vec = (pg_per_vec+1) * 4096;

  pixel_type* pixel_data = 0;
  posix_memalign((void**)&pixel_data, 4096, num_random_vector * sz_per_vec);
  assert(pixel_data != 0);
  memset((void*)pixel_data, 0, num_random_vector * sz_per_vec);

  vector_xy<pixel_type>* work = 0;
  posix_memalign((void**)&work, 4096, batch_sz * sizeof(vector_xy<pixel_type>));
  assert(work != 0);
  memset((void*)work, 0, batch_sz * sizeof(vector_xy<pixel_type>));

  // pre-generate vectors
  std::cout << "Generate " << batch_sz << std::endl;
  const auto gen_start = utility::elapsed_time_sec();
  for (int j = 0; j < batch_sz; j++) {
    const double z_intercept = z_start_dist(rnd_engine);
    const uint64_t x_intercept = x_start_dist(rnd_engine);
    const uint64_t y_intercept = y_start_dist(rnd_engine);
    const double x_slope = x_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);
    const double y_slope = y_beta_dist(rnd_engine) * 2 * (plus_or_minus(rnd_engine) ? -1 : 1);
    work[j].set(x_slope, x_intercept, y_slope, y_intercept, z_intercept, &pixel_data[j * (sz_per_vec / sizeof(pixel_type))]);
  }
  total_execution_time += utility::elapsed_time_sec(gen_start);
  std::cout << "Gen time:" << total_execution_time << std::endl;

  // fetch pixel data for each vector
  // note: this does ALL vectors in a single pass. this method is constrained
  //       by the amount of system memory available.
  std::cout << "PF " << batch_sz << std::endl;
  const auto pf_start = utility::elapsed_time_sec();
  for (uint64_t k = 0; k < size_k; k++) {
    const double time_offset = cube.timestamp(k) - cube.timestamp(0);
    #pragma omp parallel for schedule(static, NPF)
    for (uint64_t j = 0; j < batch_sz; j++) {
      // the work performed by this loop was out of cube_iterator_with_vector in vector.hpp
      vector_xy<pixel_type>& v = work[j];
      const auto xy = v.position(time_offset);
      if (time_offset >= v.z_intercept && !cube.out_of_range(xy.first, xy.second, k)) {
        v.pixels[v.npixels++] = *cube.get_pixel_addr(xy.first, xy.second, k);
      }
    }
    if (k % 10 == 0) { std::cout << k << " / " << size_k << std::endl; }
  }

  double pf_time = utility::elapsed_time_sec(pf_start);
  total_execution_time += pf_time;
  std::cout << "PF time:" << pf_time << std::endl;

  // do the median calculation
  std::cout << "Calc " << batch_sz << std::endl;
  const auto calc_start = utility::elapsed_time_sec();
  #pragma omp parallel for schedule(static, NPF)
  for (int j = 0; j < batch_sz; j++) {
    pixel_type* begin = work[j].pixels;
    pixel_type* end = &(work[j].pixels[work[j].npixels]);
    uint64_t size = (end - begin) / sizeof(pixel_type);
    // calculate median of this vector's pixel data
    std::sort(begin, end);
    uint64_t mid = size / 2;
    if (size % 2) {
      result[j].first = begin[mid];
    } else {
      pixel_type med = begin[mid];
      med += begin[mid+1];
      result[j].first = med / 2;
    }
    result[j].second = (work[j]);
  }
  double calc_time = utility::elapsed_time_sec(calc_start);
  total_execution_time += calc_time;
  std::cout << "Calc time:" << calc_time << std::endl;

  free(pixel_data);
  free(work);
  return std::make_pair(total_execution_time, result);
}

void print_top_median(const cube<pixel_type> &cube,
                      const size_t num_top,
                      std::vector<std::pair<pixel_type, vector_xy<pixel_type>>> &result) {

  // Sort the results by the descending order of median value
  std::sort(result.begin(), result.end(),
            [](const std::pair<pixel_type, vector_xy<pixel_type>> &lhd,
               const std::pair<pixel_type, vector_xy<pixel_type>> &rhd) {
              return (lhd.first > rhd.first);
            });

  // Print out the top 'num_top' median values and corresponding pixel values
  std::cout << "Top " << num_top << " median and pixel values (skip NaN value)" << std::endl;
  for (size_t i = 0; i < num_top; ++i) {
    const pixel_type median = result[i].first;
    const vector_xy<pixel_type> vector = result[i].second;

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
  uint64_t n_vec = (num_random_vector < NGEN) ? NGEN : num_random_vector;

  const auto start = utility::elapsed_time_sec();
  auto result = shoot_vector(cube, n_vec);
  double txt = utility::elapsed_time_sec(start);
  double thread_exec = result.first;

  std::cout << "#of vectors = " << n_vec
            << "\nexecution time (sec) = " << txt
            << "\nvectors/sec = " << static_cast<double>(n_vec) / txt << std::endl;

//   print_top_median(cube, std::min(num_random_vector, static_cast<size_t>(10)), result.second);

  utility::umap_fits_file::PerFits_free_cube(image_data);

  return 0;
}
