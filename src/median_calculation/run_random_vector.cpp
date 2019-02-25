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
#include <math.h>

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
#include "custom_distribution.hpp"

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

std::vector<double> read_exposuretime(const size_t size_k) {
  std::vector<double> exposuretime_list;

  const char *exposuretime_file_name = std::getenv("EXPOSURETIME_FILE");
  if (exposuretime_file_name != nullptr) {
    std::ifstream ifs(exposuretime_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open " << exposuretime_file_name << std::endl;
      std::abort();
    }
    for (double exposuretime; ifs >> exposuretime;) {
      exposuretime_list.emplace_back(exposuretime);
    }
    if (exposuretime_list.size() != size_k) {
      std::cerr << "#of lines in " << exposuretime_file_name << " is not the same as #of fits files" << std::endl;
      std::abort();
    }
  } else {
    // If a list of exposure times is not given, assume that each exposure time is 40 s
    exposuretime_list.resize(size_k);
    for (size_t i = 0; i < size_k; ++i) exposuretime_list[i] = 40;
  }

  return exposuretime_list;
}

std::vector<double> read_psf(const size_t size_k) {
  std::vector<double> psf_list;

  const char *psf_file_name = std::getenv("PSF_FILE");
  if (psf_file_name != nullptr) {
    std::ifstream ifs(psf_file_name);
    if (!ifs.is_open()) {
      std::cerr << "Cannot open " << psf_file_name << std::endl;
      std::abort();
    }
    for (double psf; ifs >> psf;) {
      psf_list.emplace_back(psf);
    }
    if (psf_list.size() != size_k) {
      std::cerr << "#of lines in " << psf_file_name << " is not the same as #of fits files" << std::endl;
      std::abort();
    }
  } else {
    // If a list of psfs is not given, assume that each psf is 1
    psf_list.resize(size_k);
    for (size_t i = 0; i < size_k; ++i) psf_list[i] = 1.0;
  }

  return psf_list;
}



	
// Function to sum across the vector, should skip nan values
template <typename iterator_type>
typename iterator_type::value_type
	vector_sum(iterator_type iterator_begin, iterator_type iterator_end) {
	using value_type = typename iterator_type::value_type;

	if (iterator_begin == iterator_end)
		return 0;

	pixel_type total = 0;
	for (auto iterator(iterator_begin); iterator != iterator_end; ++iterator) {
		const value_type value = *iterator;
		if (is_nan(value)) continue;
		total += value;
	}

	return total;
}



std::pair<double, std::vector<std::vector<pixel_type, pixel_type, vector_xy>>>
shoot_vector(const cube<pixel_type> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::vector<pixel_type, pixel_type, vector_xy>> result(num_random_vector);

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
    std::uniform_int_distribution<int> x_start_dist(0, std::get<0>(cube.size()) - 1);
    std::uniform_int_distribution<int> y_start_dist(0, std::get<1>(cube.size()) - 1);
    const char *slope_filename = std::getenv("SLOPE_PDF_FILE");
    custom_distribution slope_distribution(slope_filename);
    std::uniform_int_distribution<int> plus_or_minus(0, 1);

    // Shoot random vectors using multiple threads
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < num_random_vector; ++i) {
      
      const std::vector<double> slopes = slope_distribution();
      const double x_slope = slopes[0];
      const double y_slope = slopes[1];
      
      const double x_intercept = x_start_dist(rnd_engine);
      const double y_intercept = y_start_dist(rnd_engine);
      
      vector_xy current_vector{x_slope, x_intercept, y_slope, y_intercept};

      cube_iterator_with_vector<pixel_type> begin(cube, current_vector, 0.0);
      cube_iterator_with_vector<pixel_type> end(cube, current_vector);

      // median calculation using Torben algorithm
      const auto start = utility::elapsed_time_sec();
      result[i][0] = torben(begin, end);
	  result[i][1] = vector_sum(begin, end);
      total_execution_time += utility::elapsed_time_sec(start);
      result[i][2] = current_vector;
    }
  }

  return std::make_pair(total_execution_time, result);
}

void print_top_median(const cube<pixel_type> &cube,
                      const size_t num_top,
                      std::vector<std::vector<pixel_type, pixel_type, vector_xy>> &result) {

  // Sort the results by the descending order of median value
  std::sort(result.begin(), result.end(),
            [](const std::vector<pixel_type, pixel_type, vector_xy> &lhd,
               const std::vector<pixel_type, pixel_type, vector_xy> &rhd) {
              return (lhd[0] > rhd[0]);
            });

  // Print out the top 'num_top' median values and corresponding pixel values
  std::cout << "Top " << num_top << " median and pixel values (skip NaN value)" << std::endl;
  for (size_t i = 0; i < num_top; ++i) {
    const pixel_type median = result[i][0];
    const vector_xy vector = result[i][2];

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

// Function to write results to a csv file in the form:
// ID | MEDIAN | SUM | X_START | Y_START | X_SLOPE | Y_SLOPE
void write_tocsv(std::vector<std::vector<pixel_type, pixel_type, vector_xy>> &result) {
	std::ofstream out("vector_output.csv");

	int id = 0;
	for (auto& row : result) {
		
		out << id << ',';
		out << row[0] << ',';
		out << row[1] << ',';
		out << row[2].x_intercept << ',';
		out << row[2].y_intercept << ',';
		out << row[2].x_slope << ',';
		out << row[2].y_slope << ',';
		out << '\n';
		++id;
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
  cube<pixel_type> cube(size_x, size_y, size_k, image_data, read_timestamp(size_k), read_exposuretime(size_k), read_psf(size_k));

  const std::size_t num_random_vector = get_num_vectors();

  std::pair<double, std::vector<std::vector<pixel_type, pixel_type, vector_xy>>> result = shoot_vector(cube, num_random_vector);

  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << result.first
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / result.first << std::endl;

  print_top_median(cube, std::min(num_random_vector, static_cast<size_t>(10)), result.second);

  write_tocsv(result.second);

  utility::umap_fits_file::PerFits_free_cube(image_data);

  return 0;
}
