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
#include "velocity_distribution.hpp"

using namespace median;

using pixel_type = float;

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
	}
	else {
		// If a list of exposure times is not given, assume that each exposure time is 40 s
		exposuretime_list.resize(size_k);
		for (size_t i = 0; i < size_k; ++i) exposuretime_list[i] = 40;
	}

	return exposuretime_list;
}

vector_xy make_vector() {
	double x_int, y_int, x_slo, y_slo;

	const char *x_intercept = std::getenv("X_INTERCEPT");
	const char *y_intercept = std::getenv("Y_INTERCEPT");
	const char *x_slope = std::getenv("X_SLOPE");
	const char *y_slope = std::getenv("Y_SLOPE");

	x_int = atof(x_intercept);
	y_int = atof(y_intercept);
	x_slo = atof(x_slope);
	y_slo = atof(y_slope);
	vector_xy vector{ x_slo,x_int,y_slo,y_int };
	return vector;
}

// Function to read data info from a csv file
// Reads timestamp, psf fwhm, (ra/dec), and background sky noise 
std::tuple<std::vector<unsigned long>, std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> read_list_csv(const size_t size_k) {
	std::vector<double> psf_list;
	std::vector<unsigned long> timestamp_list;
	std::vector<std::vector<double>> ra_dec_list;
	std::vector<double> noise_list;

	const char *list_file_name = std::getenv("DATA_LIST_FILE");
	if (list_file_name != nullptr) {
		std::ifstream ifs(list_file_name);

		if (!ifs.is_open()) {
			std::cerr << "Cannot open " << list_file_name << std::endl;
			std::abort();
		}

		std::string line;
		std::getline(ifs, line); //skip first row of header info

		int check = 0;
		double psf, mjd, ra, dec, noise;
		char vala, valb, valc, vald, vale;
		std::string frame;
		double mjd_start = 0;

		while (std::getline(ifs, line))
		{
			std::istringstream iss{ line };

			std::getline(iss, frame, ',');
			do {

				if (iss >> psf >> valb >> mjd >> valc >> ra >> vald >> dec >> vale >> noise)
				{

					if (check == 0)
						mjd_start = mjd;
					unsigned long time = std::round((mjd - mjd_start) * 24 * 60 * 60 * 100); // Hundreths of a second
					std::vector<double> ra_dec = { ra, dec };

					if (psf == 0) // For nan rows we set values to the previous row
					{
						time = timestamp_list[(check - 1)];
						psf = psf_list[(check - 1)];
						ra_dec = ra_dec_list[(check - 1)];
						noise = noise_list[(check - 1)];
					}

					timestamp_list.push_back(time);
					psf_list.push_back(psf);
					ra_dec_list.push_back(ra_dec);
					noise_list.push_back(noise);
					++check;
				}
			} while (!iss.eof());
		}

		if (psf_list.size() != size_k) {
			std::cerr << "#of lines in " << list_file_name << " is not the same as #of fits files" << std::endl;
			std::abort();
		}
	}
	else {
		timestamp_list.resize(size_k);
		psf_list.resize(size_k);
		ra_dec_list.resize(size_k);
		noise_list.resize(size_k);
		for (size_t i = 0; i < size_k; ++i) {
			timestamp_list[i] = (double)i*100.0;
			psf_list[i] = 2.0;
			ra_dec_list[i] = { 0.0,0.0 };
			noise_list[i] = 100.0;
		}
	}

	return std::make_tuple(timestamp_list, psf_list, ra_dec_list, noise_list);
}

template <typename iterator_type>
std::vector<std::tuple<size_t, ssize_T, ssize_t, int, value_type>> vector_pos_info(iterator_type iterator_begin, iterator_type iterator_end) {
	using value_type = typename iterator_type::value_type;

	std::vector<std::tuple<size_t, ssize_T, ssize_t, int, value_type>> vector_info;

	for (auto iterator(iterator_begin); iterator != iterator_end; ++iterator) {
		std::tuple<pixel_type, int, ssize_t, ssize_t, size_t> vec_info = iterator.vector_pos_info();
		const value_type value = std::get<0>(vec_info);
		int num_pixels = std::get<1>(vec_info);
		ssize_t x_pos = std::get<2>(vec_info);
		ssize_t y_pos = std::get<3>(vec_info);
		size_t k_pos = std::get<4>(vec_info);
		auto tup = std::make_tuple(k_pos, x_pos, y_pos, num_pixels, value);

		vector_info.push_back(tup);		
	}



	return vector_info;
}

// Function to write results to a csv file in the form:
// K position | X position | Y position | Number of pixels hit | SUM 
void write_tocsv(std::vector<std::tuple<size_t, ssize_T, ssize_t, int, value_type>> &result) {
	std::ofstream out("vector_check_output.csv");

	out << "K_POS,X_POS,Y_POS,NUM_PIXELS,SUM\n";

	for (auto& row : result) {

		out << std::get<0>(row) << ',';
		out << std::get<1>(row) << ',';
		out << std::get<2>(row) << ',';
		out << std::get<3>(row) << ',';
		out << std::get<4>(row) << ',';
		out << '\n';
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
	std::tuple<std::vector<unsigned long>, std::vector<double>, std::vector<std::vector<double>>, std::vector<double>> lists = read_list_csv(size_k);
	cube<pixel_type> cube(size_x, size_y, size_k, image_data, std::get<0>(lists), read_exposuretime(size_k), std::get<1>(lists), std::get<2>(lists), std::get<3>(lists));

	//get intercepts and slopes for a given vector
	vector_xy current_vector = make_vector();

	cube_iterator_with_vector<pixel_type> begin(cube, current_vector, 0.0);
	cube_iterator_with_vector<pixel_type> end(cube, current_vector);

	auto vector_info = vector_pos_info(begin, end);

	write_tocsv(vector_info);

	//save image cutout of each frame?


	return 0;
}