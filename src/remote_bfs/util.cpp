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
#include "util.hpp"

#include "umap/umap.h"
#include "../utility/mmap.hpp"

void usage() {
  std::cout << "BFS options:"
            << "-n\t#vertices\n"
            << "-m\t#edges\n"
            << "-g\tGraph file name\n"
            << "-s\tUse system mmap" << std::endl;
}

void parse_options(int argc, char **argv,
                   bfs_options &options) {
  int c;
  while ((c = getopt(argc, argv, "n:m:g:sh")) != -1) {
    switch (c) {
      case 'n': /// Required
        options.num_vertices = std::stoull(optarg);
        break;

      case 'm': /// Required
        options.num_edges = std::stoull(optarg);
        break;

      case 'g': /// Required
        options.graph_file_name = optarg;
        break;

      case 's':options.is_server = true;
        break;

      case 'h':usage();
        std::exit(0);
    }
  }
}


void disp_umap_env_variables() {
  std::cout
      << "Environment Variable Configuration (command line arguments obsolete):\n"
      << "UMAP_PAGESIZE                   - currently: " << umapcfg_get_umap_page_size() << " bytes\n"
      << "UMAP_PAGE_FILLERS               - currently: " << umapcfg_get_num_fillers() << " fillers\n"
      << "UMAP_PAGE_EVICTORS              - currently: " << umapcfg_get_num_evictors() << " evictors\n"
      << "UMAP_READ_AHEAD                 - currently: " << umapcfg_get_read_ahead() << " pages\n"
      << "UMAP_BUFSIZE                    - currently: " << umapcfg_get_max_pages_in_buffer() << " pages\n"
      << "UMAP_EVICT_LOW_WATER_THRESHOLD  - currently: " << umapcfg_get_evict_low_water_threshold() << " percent full\n"
      << "UMAP_EVICT_HIGH_WATER_THRESHOLD - currently: " << umapcfg_get_evict_high_water_threshold()
      << " percent full\n"
      << std::endl;
}


void disp_bfs_options(const bfs_options &options) {
  std::cout << "BFS options:"
            << "\n#vertices: " << options.num_vertices
            << "\n#edges: " << options.num_edges
            << "\nGraph file: " << options.graph_file_name
            << "\nis Server: " << options.is_server << std::endl;
}

std::string disp_omp_schedule(const int kind_in_int) {
#ifdef _OPENMP
  if (kind_in_int == omp_sched_static) {
    return std::string("omp_sched_static (" + std::to_string(kind_in_int) + ")");
  } else if (kind_in_int == omp_sched_dynamic) {
    return std::string("omp_sched_dynamic (" + std::to_string(kind_in_int) + ")");
  } else if (kind_in_int == omp_sched_guided) {
    return std::string("omp_sched_guided (" + std::to_string(kind_in_int) + ")");
  } else if (kind_in_int == omp_sched_auto) {
    return std::string("omp_sched_auto (" + std::to_string(kind_in_int) + ")");
  }
  return std::string("Unknown kind (" + std::to_string(kind_in_int) + ")");
#else
  return std::string("OpenMP is not supported");
#endif
}


void print_num_page_faults() {
  const auto num_page_faults = utility::get_num_page_faults();
  std::cout << "#of minor page faults\t" << num_page_faults.first << std::endl;
  std::cout << "#of major page faults\t" << num_page_faults.second << std::endl;
}
