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
#include <tuple>
#include <fstream>

#include "util.hpp"
#include "../utility/umap_file.hpp"
#include "../utility/time.hpp"
#include "../utility/file.hpp"


std::pair<uint64_t *, uint64_t *>
map_graph(const bfs_options &options) {

  /* round up to page aligned size */
  const size_t original_size  = (options.num_vertices + 1 + options.num_edges) * sizeof(uint64_t);
  const size_t umap_page_size = umapcfg_get_umap_page_size();
  const size_t aligned_pages  = (original_size - 1)/umap_page_size + 1;
  const size_t aligned_size   = umap_page_size * aligned_pages;

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


int main(int argc, char **argv) {
  bfs_options options;

  parse_options(argc, argv, options);
  disp_bfs_options(options);
  if (!options.use_mmap) disp_umap_env_variables();

  std::cout << "Initial #of page faults" << std::endl;
  print_num_page_faults();

  const uint64_t *index = nullptr;
  const uint64_t *edges = nullptr;
  std::tie(index, edges) = map_graph(options);

  // Array to store each vertex's level (a distance from the source vertex)
  std::vector<uint16_t> level(options.num_vertices);

  // bitmap data to store 'visited' information
  std::vector<uint64_t> visited_filter(utility::bitmap_size(options.num_vertices));

  init_bfs(options.num_vertices, level.data(), visited_filter.data());
  find_bfs_root(options.num_vertices, index, level.data());

  std::cout << "Before BFS #of page faults" << std::endl;
  print_num_page_faults();
  const auto bfs_start_time = utility::elapsed_time_sec();
  const uint16_t max_level = run_bfs(options.num_vertices, index, edges, level.data(), visited_filter.data());
  const auto bfs_time = utility::elapsed_time_sec(bfs_start_time);
  std::cout << "BFS took (s)\t" << bfs_time << std::endl;
  std::cout << "After BFS #of page faults" << std::endl;
  print_num_page_faults();

  count_level(options.num_vertices, max_level, level.data());

  utility::unmap_file(options.use_mmap,
                      utility::get_file_size(options.graph_file_name),
                      const_cast<uint64_t *>(index));

  return 0;
}


