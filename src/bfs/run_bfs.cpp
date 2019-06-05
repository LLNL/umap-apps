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
#include <unistd.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <fstream>

#include "umap/umap.h"
#include "bfs_kernel.hpp"
#include "../utility/umap_file.hpp"
#include "../utility/bitmap.hpp"
#include "../utility/time.hpp"
#include "../utility/file.hpp"

struct bfs_options {
  size_t num_vertices{0};
  size_t num_edges{0};
  std::string graph_file_name;
  bool use_mmap{false};
};

void disp_umap_env_variables() {
  std::cout
      << "Environment Variable Configuration (command line arguments obsolete):\n"
      << "UMAP_PAGESIZE                   - currently: " << umapcfg_get_umap_page_size() << " bytes\n"
      << "UMAP_PAGE_FILLERS               - currently: " << umapcfg_get_num_fillers() << " fillers\n"
      << "UMAP_PAGE_EVICTORS              - currently: " << umapcfg_get_num_evictors() << " evictors\n"
      << "UMAP_READ_AHEAD                 - currently: " << umapcfg_get_read_ahead() << " evictors\n"
      << "UMAP_BUFSIZE                    - currently: " << umapcfg_get_max_pages_in_buffer() << " pages\n"
      << "UMAP_EVICT_LOW_WATER_THRESHOLD  - currently: " << umapcfg_get_evict_low_water_threshold() << " percent full\n"
      << "UMAP_EVICT_HIGH_WATER_THRESHOLD - currently: " << umapcfg_get_evict_high_water_threshold()
      << " percent full\n"
      << std::endl;
}

void usage() {
  std::cout << "BFS options:"
            << "-n\t#vertices\n"
            << "-m\t#edges\n"
            << "-g\tGraph file name\n"
            << "-s\tUse system mmap" << std::endl;

  disp_umap_env_variables();
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

      case 's':options.use_mmap = true;
        break;

      case 'h':usage();
        std::exit(0);
    }
  }
}

void disp_bfs_options(const bfs_options &options) {
  std::cout << "BFS options:"
            << "\n#vertices: " << options.num_vertices
            << "\n#edges: " << options.num_edges
            << "\nGraph file: " << options.graph_file_name
            << "\nUse system mmap: " << options.use_mmap << std::endl;
}

size_t calculate_umap_pagesize_aligned_graph_file_size(const size_t num_vertices, const size_t num_edges) {
  const size_t original_size = (num_vertices + 1 + num_edges) * sizeof(uint64_t);
  const size_t umap_page_size = umapcfg_get_umap_page_size();
  const size_t aligned_graph_size = (original_size % umap_page_size == 0)
                                    ? original_size
                                    : (original_size + (umap_page_size - original_size % umap_page_size));
  return aligned_graph_size;
}

std::pair<uint64_t *, uint64_t *>
map_graph(const bfs_options &options) {

  // Umap requires a pagesize aligned file
  if (!options.use_mmap) {
    const size_t size = calculate_umap_pagesize_aligned_graph_file_size(options.num_vertices, options.num_edges);
    if (!utility::extend_file_size(options.graph_file_name, size)) {
      std::cerr << "Failed to extend the graph file to " << size << std::endl;
      std::abort();
    }
  }

  void *const map_raw_address = utility::map_in_file(options.graph_file_name,
                                                     false,
                                                     true,
                                                     options.use_mmap,
                                                     utility::get_file_size(options.graph_file_name),
                                                     nullptr);
  if (!map_raw_address) {
    std::cerr << "Failed to map the graph" << std::endl;
    std::abort();
  }

  uint64_t *const index = static_cast<uint64_t *>(map_raw_address);
  const uint64_t edges_offset = options.num_vertices + 1;
  uint64_t *const edges = static_cast<uint64_t *>(map_raw_address) + edges_offset;

  return std::make_pair(index, edges);
}

void find_bfs_root(const size_t num_vertices, const uint64_t *const index, uint16_t *const level) {
  for (uint64_t src = 0; src < num_vertices; ++src) {
    const size_t degree = index[src + 1] - index[src];
    if (degree > 0) {
      level[src] = 0;
      std::cout << "BFS root: " << src << std::endl;
      return;
    }
  }
  std::cerr << "Can not find a proper root vertex; all vertices do not have any edges?" << std::endl;
  std::abort();
}

void count_level(const size_t num_vertices, const uint16_t max_level, const uint16_t *const level) {

  std::vector<size_t> cnt(max_level + 1, 0);
  for (uint64_t i = 0; i < num_vertices; ++i) {
    if (level[i] == bfs::k_infinite_level) continue;
    if (level[i] > max_level) {
      std::cerr << "Invalid level: " << level[i] << " > " << max_level << std::endl;
      return;
    }
    ++cnt[level[i]];
  }

  std::cout << "Level\t#vertices" << std::endl;
  for (uint16_t i = 0; i <= max_level; ++i) {
    std::cout << i << "\t" << cnt[i] << std::endl;
  }
}

int main(int argc, char **argv) {
  bfs_options options;

  parse_options(argc, argv, options);
  disp_bfs_options(options);
  if (!options.use_mmap) disp_umap_env_variables();

  const uint64_t *index = nullptr;
  const uint64_t *edges = nullptr;
  std::tie(index, edges) = map_graph(options);

  // Array to store each vertex's level (a distance from the source vertex)
  std::vector<uint16_t> level(options.num_vertices);

  // bitmap data to store 'visited' information
  std::vector<uint64_t> visited_filter(utility::bitmap_size(options.num_vertices));

  bfs::init_bfs(options.num_vertices, level.data(), visited_filter.data());
  find_bfs_root(options.num_vertices, index, level.data());

  const auto bfs_start_time = utility::elapsed_time_sec();
  const uint16_t max_level = bfs::run_bfs(options.num_vertices, index, edges, level.data(), visited_filter.data());
  const auto bfs_time = utility::elapsed_time_sec(bfs_start_time);
  std::cout << "BFS took (s) " << bfs_time << std::endl;

  count_level(options.num_vertices, max_level, level.data());

  utility::unmap_file(options.use_mmap,
                      utility::get_file_size(options.graph_file_name),
                      const_cast<uint64_t *>(index));

  return 0;
}
