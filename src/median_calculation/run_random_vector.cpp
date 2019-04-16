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
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utility/commandline.hpp"
#include "../utility/umap_fits_file.hpp"
#include "../utility/time.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include "velocity_distribution.hpp"

using namespace median;
using pixel_type = float;
using fits_cube = utility::umap_fits_file::umap_fits_cube<pixel_type>;
constexpr size_t default_num_random_vector = 100000;

std::size_t get_num_vectors() {
  std::size_t num_random_vector = default_num_random_vector;
  const char *buf = std::getenv("NUM_VECTORS");
  if (buf != nullptr) {
    num_random_vector = std::stoll(buf);
  }
  return num_random_vector;
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

// Function to calculate relevant information about a given vector
// Returns: <SNR, weighted sum, number of frames intersected>
template <typename iterator_type>
std::tuple<double, typename iterator_type::value_type, int> vector_info(iterator_type iterator_begin, iterator_type iterator_end) {
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
    std::tuple<pixel_type, int, double, double> snr_info = iterator.snr_info();
    const value_type value = std::get<0>(snr_info);
    int num_pixels = std::get<1>(snr_info);

    if (num_pixels == 0) continue; // For when the vector hits nothing in an image

    total_signal += value;
    ++frame_num;

    // SNR calculation
    double B = std::get<3>(snr_info); // pull background noise from list
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

std::vector<std::tuple<vector_xy, double, double, int>> shoot_vector(const fits_cube &d_cube, const std::size_t num_random_vector) {
  // Array to store results of the median calculation
  std::vector<std::tuple<vector_xy, double, double, int>> result(num_random_vector);
  uint64_t rng_time = 0;
  uint64_t vector_time = 0;
  uint64_t total_time = 0;
  int num_threads = 1;

  const auto start = utility::elapsed_time();
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
#ifdef _OPENMP
    if (omp_get_thread_num() == 0) {
      num_threads == omp_get_num_threads();
    }
    std::mt19937 rnd_engine(123 + omp_get_thread_num());
#else
    std::mt19937 rnd_engine(123);
#endif
//     std::uniform_int_distribution<std::size_t> start_dist(0, d_cube.cube_size() - 1);

    // Generate a slope distribution from a given file or a beta distribution if no file given
    const char *slope_filename = std::getenv("SLOPE_PDF_FILE");

    // Function takes both potential filename, (optional) a/b arguments for beta distribution, and pixel scale
    slope_distribution slope_dist(slope_filename,3,2,0.26);

    // Shoot random vectors using multiple threads
#ifdef _OPENMP
#pragma omp for
#endif
    for (int i = 0; i < num_random_vector; ++i) {
      double tmp = 0;

      // record how much time was spent generating random vectors
      const auto rng_start = utility::elapsed_time();
        std::vector<double> slopes = slope_dist(rnd_engine, d_cube.ra_dec(0)[0], d_cube.ra_dec(0)[1]);
        double x_slope = slopes[0];
        double y_slope = slopes[1];

        auto intercept = d_cube.get_rnd_coord(rnd_engine, x_slope, y_slope);
        double x_intercept = std::get<0>(intercept);
        double y_intercept = std::get<1>(intercept);
      __sync_fetch_and_add(&rng_time, utility::elapsed_time_nsec(rng_start));

      // record how much time was spent tracing the vectors
      const auto vec_start = utility::elapsed_time();
        vector_xy current_vector{x_slope, x_intercept, y_slope, y_intercept};
        cube_iterator_with_vector<pixel_type> begin(d_cube, current_vector, 0.0);
        cube_iterator_with_vector<pixel_type> end(d_cube, current_vector);

      // vector info stored as [VECTOR_XY, SNR, SUM, NUMBER OF FRAMES]
        std::tuple<double, double, int> v_info = vector_info(begin, end);
        result[i] = std::make_tuple(current_vector,std::get<0>(v_info),std::get<1>(v_info),std::get<2>(v_info));
      __sync_fetch_and_add(&vector_time, utility::elapsed_time_nsec(vec_start));
    }
  }

  total_time += utility::elapsed_time_nsec(start);
  std::cout << "total: " << ((double)total_time) / 1e9 << " rng: " << ((double)rng_time) / 1e9 << " vec: " << ((double)vector_time) / 1e9 << std::endl;
  return result;
}

// Function to write results to a csv file in the form:
// ID | X_INTERCEPT | Y_INTERCEPT | X_SLOPE | Y_SLOPE | SNR | SUM | NUMBER OF FRAMES HIT
void write_tocsv(std::vector<std::tuple<vector_xy, double, double, int>> &result) {
  std::ofstream out("vector_output.csv");
  out << "ID,X_INTERCEPT,Y_INTERCEPT,X_SLOPE,Y_SLOPE,SNR,SUM,NUMBER_OF_FRAMES_HIT\n";
  long long id = 1;
  for (auto& row : result) {
    out << id << ','
        << std::get<0>(row).x_intercept << ','
        << std::get<0>(row).y_intercept << ','
        << std::get<0>(row).x_slope << ','
        << std::get<0>(row).y_slope << ','
        << std::get<1>(row) << ','
        << std::get<2>(row) << ','
        << std::get<3>(row) << std::endl;
    ++id;
  }
}

int main(int argc, char **argv) {
  utility::umt_optstruct_t options;
  umt_getoptions(&options, argc, argv);

#ifdef _OPENMP
  omp_set_num_threads(options.numthreads);
#endif

  const char *data_list_name = std::getenv("DATA_LIST_FILE");
  const char *exptime_list_name = std::getenv("EXPOSURETIME_FILE");
  fits_cube d_cube(options.dirname, data_list_name, exptime_list_name);
  const std::size_t num_random_vector = get_num_vectors();

  const auto main_start = utility::elapsed_time_sec();
  auto result = shoot_vector(d_cube, num_random_vector);
  const double main_exec_time = utility::elapsed_time_sec(main_start);

  std::cout << "#of vectors = " << num_random_vector
            << "\nexecution time (sec) = " << main_exec_time
            << "\nvectors/sec = " << static_cast<double>(num_random_vector) / main_exec_time << std::endl;

  write_tocsv(result);
  return 0;
}
