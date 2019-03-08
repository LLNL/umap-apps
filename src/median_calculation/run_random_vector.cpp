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
#include "utility.hpp"
#include "vector.hpp"
#include "cube.hpp"
#include "beta_distribution.hpp"
#include "custom_distribution.hpp"
#include "distribution_test.hpp"

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
    for (size_t i = 0; i < size_k; ++i) psf_list[i] = 2.5;
  }

  return psf_list;
}



	
// Function to calculate relevant information about a given vector
// Returns: <SNR, weighted sum, number of frames intersected>
template <typename iterator_type>
std::tuple<double, typename iterator_type::value_type, int>
	vector_info(iterator_type iterator_begin, iterator_type iterator_end) {
	using value_type = typename iterator_type::value_type;

	if (iterator_begin == iterator_end)
		return std::tuple<double, value_type, int> (0,0,0);

	// DECAM info
	double dark_noise = 0.417; // electrons per pixel per second
	double readout_noise = 7; // electrons

	value_type total_signal = 0;
	double total_B = 0;
	double total_R = 0;
	double total_D = 0;
	int frame_num = 0;

	for (auto iterator(iterator_begin); iterator != iterator_end; ++iterator) {
		std::tuple<pixel_type, int, double> snr_info = iterator.snr_info();
		const value_type value = std::get<0>(snr_info);
		int num_pixels = std::get<1>(snr_info);
		
		if (num_pixels == 0) continue; // For when the vector hits nothing in an image

		total_signal += value;
		++frame_num;

		// SNR calculation
		double B = 10*dark_noise; // pull background noise from list???
		double exp_time = std::get<2>(snr_info);

		total_B += B * num_pixels * exp_time;
		total_R += num_pixels * readout_noise*readout_noise;
		total_D += dark_noise * num_pixels * exp_time;
	}	

	double SNR = 0;
	if ((total_signal > 0) && (frame_num != 0))
		SNR = total_signal/sqrt(total_signal + total_B + total_D + total_R);
		
	return std::tuple<double, value_type, int> (SNR, total_signal, frame_num);
}


std::pair<double, std::vector<std::tuple<vector_xy, double, double, int>>>
shoot_vector(const cube<pixel_type> &cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::tuple<vector_xy, double, double, int>> result(num_random_vector);

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

	// Generate a slope distribution from a given file or a beta distribution if no file given
    
    const char *slope_filename = std::getenv("SLOPE_PDF_FILE");

	// Function takes both potential filename and potential (optional) a/b arguments for beta distribution
    slope_distribution slope_dist(slope_filename,3,2);

    
    // Shoot random vectors using multiple threads
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < num_random_vector; ++i) {
      
      std::vector<double> slopes = slope_dist(rnd_engine);
      double x_slope = slopes[0];
      double y_slope = slopes[1];
      double x_intercept = x_start_dist(rnd_engine);
      double y_intercept = y_start_dist(rnd_engine);
      
      vector_xy current_vector{x_slope, x_intercept, y_slope, y_intercept};

      cube_iterator_with_vector<pixel_type> begin(cube, current_vector, 0.0);
      cube_iterator_with_vector<pixel_type> end(cube, current_vector);
      	
	  // vector info stored as [VECTOR_XY, SNR, SUM, NUMBER OF FRAMES]
      const auto start = utility::elapsed_time_sec();
	  std::tuple<double, double, int> v_info = vector_info(begin, end);
	  result[i] = std::make_tuple(current_vector,std::get<0>(v_info),std::get<1>(v_info),std::get<2>(v_info));
      total_execution_time += utility::elapsed_time_sec(start);
    }
  }

  return std::make_pair(total_execution_time, result);
}

/*
void print_top_median(const cube<pixel_type> &cube,
                      const size_t num_top,
                      std::vector<std::tuple<pixel_type, pixel_type, vector_xy>> &result) {

  // Sort the results by the descending order of median value
  std::sort(result.begin(), result.end(),
            [](const std::tuple<pixel_type, pixel_type, vector_xy> &lhd,
               const std::tuple<pixel_type, pixel_type, vector_xy> &rhd) {
              return (std::get<0>(lhd) > std::get<0>(rhd));
            });

  // Print out the top 'num_top' median values and corresponding pixel values
  std::cout << "Top " << num_top << " median and pixel values (skip NaN value)" << std::endl;
  for (size_t i = 0; i < num_top; ++i) {
    const pixel_type median = std::get<0>(result[i]);
    const vector_xy vector = std::get<2>(result[i]);

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
*/

// Function to write results to a csv file in the form:
// ID | X_INTERCEPT | Y_INTERCEPT | X_SLOPE | Y_SLOPE | SNR | SUM | NUMBER OF FRAMES HIT
void write_tocsv(std::vector<std::tuple<vector_xy, double, double, int>> &result) {
	std::ofstream out("vector_output.csv");

	out << "ID,X_INTERCEPT,Y_INTERCEPT,X_SLOPE,Y_SLOPE,SNR,SUM,NUMBER_OF_FRAMES_HIT\n";

	long long id = 1;
	for (auto& row : result) {
		
		out << id << ',';
		out << std::get<0>(row).x_intercept << ',';
		out << std::get<0>(row).y_intercept << ',';
		out << std::get<0>(row).x_slope << ',';
		out << std::get<0>(row).y_slope << ',';
		out << std::get<1>(row) << ',';
		out << std::get<2>(row) << ',';
		out << std::get<3>(row) << ',';
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

  auto result = shoot_vector(cube, num_random_vector);
  
  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << result.first
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / result.first << std::endl;

  //print_top_median(cube, std::min(num_random_vector, static_cast<size_t>(10)), result.second);

  write_tocsv(result.second);

  utility::umap_fits_file::PerFits_free_cube(image_data);

  return 0;
}
