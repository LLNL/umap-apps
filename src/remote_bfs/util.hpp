#ifndef BFS_UTIL_H
#define BFS_UTIL_H

#include <cstddef>
#include <vector>
#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utility/bitmap.hpp"


struct bfs_options {
  size_t num_vertices{0};
  size_t num_edges{0};
  std::string graph_file_name;
  bool use_mmap{false};
};

void usage();

void parse_options(int argc, char **argv,
                   bfs_options &options);

void disp_umap_env_variables();

void disp_bfs_options(const bfs_options &options);
  
std::string disp_omp_schedule(const int kind_in_int);

void print_num_page_faults();

void print_num_page_faults();
  
void init_bfs(const size_t num_vertices,
	      uint16_t *const level,
	      uint64_t *visited_filter);

uint16_t run_bfs(const size_t num_vertices,
                 const uint64_t *const index,
                 const uint64_t *const edges,
                 uint16_t *const level,
                 uint64_t *visited_filter);

void find_bfs_root(const size_t num_vertices,
		   const uint64_t *const index,
		   uint16_t *const level);

void count_level(const size_t num_vertices,
		 const uint16_t max_level,
		 const uint16_t *const level);

#endif
