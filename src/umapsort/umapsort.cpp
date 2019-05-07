//////////////////////////////////////////////////////////////////////////////
// Copyright 2017-2019 Lawrence Livermore National Security, LLC and other
// UMAP Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: LGPL-2.1-only
//////////////////////////////////////////////////////////////////////////////
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif // _GNU_SOURCE

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <random>
#include <string>
#include <vector>
#include <parallel/algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <omp.h>

#include "umap/umap.h"
#include "../utility/commandline.hpp"
#include "../utility/umap_file.hpp"
#include "../utility/time.hpp"

using namespace std;

bool sort_ascending = true;

void initdata(uint64_t *region, uint64_t rlen) {
  fprintf(stderr, "initdata: %p, from %lu to %lu\n", region, (rlen), (rlen - rlen));
#pragma omp parallel for
  for(uint64_t i=0; i < rlen; ++i)
    region[i] = (uint64_t) (rlen - i);
}

void validatedata(uint64_t *region, uint64_t rlen) {
  if (sort_ascending == true) {
#pragma omp parallel for
    for(uint64_t i = 0; i < rlen; ++i) {
        if (region[i] != (i+1)) {
            fprintf(stderr, "Worker %d found an error at index %lu, %lu != lt %lu!\n",
                            omp_get_thread_num(), i, region[i], i+1);

            if (i < 3) {
                fprintf(stderr, "Context ");
                for (int j=0; j < 7; j++) {
                    fprintf(stderr, "%lu ", region[j]);
                }
                fprintf(stderr, "\n");
            }
            else if (i > (rlen-4)) {
                fprintf(stderr, "Context ");
                for (uint64_t j=rlen-8; j < rlen; j++) {
                    fprintf(stderr, "%lu ", region[j]);
                }
                fprintf(stderr, "\n");
            }
            else {
                fprintf(stderr,
                    "Context i-3 i-2 i-1 i i+1 i+2 i+3:%lu %lu %lu %lu %lu %lu %lu\n",
                    region[i-3], region[i-2], region[i-1], region[i], region[i+1], region[i+2], region[i+3]);
            }
        }
    }
  }
  else {
#pragma omp parallel for
    for(uint64_t i = 0; i < rlen; ++i) {
        if(region[i] != (rlen - i)) {
            fprintf(stderr, "Worker %d found an error at index %lu, %lu != %lu!\n",
                            omp_get_thread_num(), i, region[i], (rlen - i));

            if (i < 3) {
                fprintf(stderr, "Context ");
                for (int j=0; j < 7; j++) {
                    fprintf(stderr, "%lu ", region[j]);
                }
                fprintf(stderr, "\n");
            }
            else if (i > (rlen-4)) {
                fprintf(stderr, "Context ");
                for (uint64_t j=rlen-8; j < rlen; j++) {
                    fprintf(stderr, "%lu ", region[j]);
                }
                fprintf(stderr, "\n");
            }
            else {
                fprintf(stderr,
                    "Context i-3 i-2 i-1 i i+1 i+2 i+3:%lu %lu %lu %lu %lu %lu %lu\n",
                    region[i-3], region[i-2], region[i-1], region[i], region[i+1], region[i+2], region[i+3]);
            }
            exit(1);
        }
    }
  }
}

int main(int argc, char **argv)
{
  utility::umt_optstruct_t options;
  uint64_t pagesize;
  uint64_t totalbytes;
  uint64_t arraysize;
  std::vector<void*> mappings;
  std::vector<std::string> filenames;
  std::vector<uint64_t> mapsize;
  void* range;

  auto start = utility::elapsed_time_sec();

  umt_getoptions(&options, argc, argv);

  pagesize = (uint64_t)utility::get_umap_page_size();

  omp_set_num_threads(options.numthreads);

  totalbytes = options.numpages*pagesize;
  range = utility::map_in_file(options.filename, options.initonly, options.noinit, options.usemmap, totalbytes);
  if (range == nullptr)
    return -1;

  if (options.numfiles > 1) {
    std::cout << "Setting up to operate on " << options.numfiles << " files." << std::endl;

    uint64_t pages_per_file = options.numpages / options.numfiles;
    uint64_t pages_in_last_file = pages_per_file + (options.numpages % options.numfiles);

    for ( int i = 0; i < options.numfiles; ++i) {
      stringstream ss;

      uint64_t pagecount = (i == (options.numfiles-1)) ? pages_per_file : pages_in_last_file;
      uint64_t bcount = pagesize * pagecount;
      mapsize.push_back(bcount);

      ss << options.filename << "." << i;
      filenames.push_back(ss.str());

      std::cout << "Mapping " << mapsize[i] << " bytes at " << range << " to " << filenames[i] << std::endl;

      void* val = utility::map_in_file(filenames[i], options.initonly, options.noinit, options.usemmap, mapsize[i], range);
      if (val == nullptr) {
        std::cerr << "Failed to map " << filenames[i] << std::endl;
        return -1;
      }
      mappings.push_back(val);
      range = (void*)((char*)range + mapsize[i]);
    }

  }
  else {
    mappings.push_back(range);
    mapsize.push_back(totalbytes);
  }

  fprintf(stderr, "umap INIT took %f seconds\n", utility::elapsed_time_sec(start));
  fprintf(stderr, "%lu pages, %lu bytes, %lu threads\n", options.numpages, totalbytes, options.numthreads);

  uint64_t *arr = (uint64_t *) mappings[0];
  arraysize = totalbytes/sizeof(uint64_t);

  start = utility::elapsed_time_sec();
  if ( !options.noinit ) {
    // init data
    initdata(arr, arraysize);
    fprintf(stderr, "INIT took %f seconds\n", utility::elapsed_time_sec(start));
  }

  if ( !options.initonly )
  {
    start = utility::elapsed_time_sec();
    sort_ascending = (arr[0] != 1);

    if (sort_ascending == true) {
      printf("Sorting in Ascending Order\n");
      __gnu_parallel::sort(arr, &arr[arraysize], std::less<uint64_t>(), __gnu_parallel::quicksort_tag());
    }
    else {
      printf("Sorting in Descending Order\n");
      __gnu_parallel::sort(arr, &arr[arraysize], std::greater<uint64_t>(), __gnu_parallel::quicksort_tag());
    }

    fprintf(stderr, "Sort took %f seconds\n", utility::elapsed_time_sec(start));

    start = utility::elapsed_time_sec();
    validatedata(arr, arraysize);
    fprintf(stderr, "Validate took %f seconds\n", utility::elapsed_time_sec(start));
  }

  start = utility::elapsed_time_sec();

  if (options.numfiles > 1) {
    for ( int i = 0; i < options.numfiles; ++i) {
      std::cout << "Unmapping " << mapsize[i] << " bytes from " << mappings[i] << std::endl;
      utility::unmap_file(options.usemmap, mapsize[i], mappings[i]);
    }
  }
  else {
    utility::unmap_file(options.usemmap, mapsize[0], mappings[0]);
  }

  fprintf(stderr, "umap TERM took %f seconds\n", utility::elapsed_time_sec(start));

  return 0;
}
