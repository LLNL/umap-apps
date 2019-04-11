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

#ifndef UMAP_TEST_LIB_UTILITY_TIME_HPP
#define UMAP_TEST_LIB_UTILITY_TIME_HPP

#include <chrono>

namespace utility {
inline std::chrono::high_resolution_clock::time_point elapsed_time() {
  return std::chrono::high_resolution_clock::now();
}

inline std::chrono::high_resolution_clock::time_point elapsed_time_sec() {
  return std::chrono::high_resolution_clock::now();
}

inline std::chrono::high_resolution_clock::time_point elapsed_time_usec() {
  return std::chrono::high_resolution_clock::now();
}

inline std::chrono::duration<unsigned long long int, std::nano> elapsed_time(const std::chrono::high_resolution_clock::time_point &tic) {
  return std::chrono::high_resolution_clock::now() - tic;
}

inline double elapsed_time_sec(const std::chrono::high_resolution_clock::time_point &tic) {
  return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(elapsed_time(tic)).count() / 1e6);
}

inline unsigned long long int elapsed_time_nsec(const std::chrono::high_resolution_clock::time_point &tic) {
  return static_cast<unsigned long long int>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed_time(tic)).count());
}
}
#endif //UMAP_TEST_LIB_UTILITY_TIME_HPP
