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
#include "../utility/mmap.hpp"

void parse_options(int argc, char **argv,
                   size_t &num_vertices, size_t &num_edges,
                   std::string &graph_file_name,
                   std::string &bfs_level_reference_file_name) {
  num_vertices = 0;
  num_edges = 0;
  graph_file_name = "";
  bfs_level_reference_file_name = "";

  int c;
  while ((c = getopt(argc, argv, "n:m:g:l:h")) != -1) {
    switch (c) {
      case 'n': /// Required
        num_vertices = std::stoull(optarg);
        break;

      case 'm': /// Required
        num_edges = std::stoull(optarg);
        break;

      case 'g': /// Required
        graph_file_name = optarg;
        break;

      case 'l': /// Required
        bfs_level_reference_file_name = optarg;
        break;

      case 'h':
        // usage();
        break;
    }
  }
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
map_graph(const size_t num_vertices, const size_t num_edges, const std::string &graph_file_name) {

  const size_t graph_size = (num_vertices + 1 + num_edges) * sizeof(uint64_t);
  /*
  int fd = -1;
  void *map_raw_address = nullptr;
  std::tie(fd, map_raw_address) = utility::map_file_read_mode(graph_file_name, nullptr, graph_size, 0);
  */
  const size_t size = calculate_umap_pagesize_aligned_graph_file_size(num_vertices, num_edges);
  if (!utility::extend_file_size(graph_file_name, size)) {
    std::cerr << "Failed to extend the graph file to " << size << std::endl;
    std::abort();
  }
  void *const map_raw_address = utility::umap_in_file(graph_file_name,
						      false, //no write, read only
						      false,
						      true,
						      false, //do NOT use mmap
						      utility::get_file_size(graph_file_name),
						      nullptr);

  if ( map_raw_address == nullptr) {
    std::cerr << "Failed to map the graph" << std::endl;
    std::abort();
  }

  uint64_t *index = static_cast<uint64_t *>(map_raw_address);
  const std::ptrdiff_t edges_offset = num_vertices + 1;
  uint64_t *edges = static_cast<uint64_t *>(map_raw_address) + edges_offset;

  return std::make_pair(index, edges);
}

void validate_level(const std::vector<uint16_t>& level, const std::string& bfs_level_reference_file_name) {

  std::ifstream ifs(bfs_level_reference_file_name);
  if (!ifs.is_open()) {
    std::cerr << "Can not open: "<< bfs_level_reference_file_name << std::endl;
    std::abort();
  }

  std::vector<uint16_t> ref_level(level.size(), bfs::k_infinite_level);
  uint64_t id;
  uint16_t lv;
  while (ifs >> id >> lv) {
    ref_level[id] = lv;
  }

  if (level != ref_level) {
    std::cerr << "BFS level is wrong" << std::endl;
    std::abort();
  }
}

int main(int argc, char **argv) {
  size_t num_vertices;
  size_t num_edges;
  std::string graph_file_name;
  std::string bfs_level_reference_file_name;

  parse_options(argc, argv, num_vertices, num_edges, graph_file_name, bfs_level_reference_file_name);

  const uint64_t *index = nullptr;
  const uint64_t *edges = nullptr;
  std::tie(index, edges) = map_graph(num_vertices, num_edges, graph_file_name);

  std::vector<uint16_t> level(num_vertices); // Array to store each vertex's level (a distance from the source vertex)
  std::vector<uint64_t> visited_filter(utility::bitmap_size(num_vertices)); // bitmap data to store 'visited' information

  bfs::init_bfs(num_vertices, level.data(), visited_filter.data());
  level[0] = 0; // Start from vertex 0

  const auto bfs_start_time = utility::elapsed_time_sec();
  bfs::run_bfs(num_vertices, index, edges, level.data(), visited_filter.data());
  const auto bfs_time = utility::elapsed_time_sec(bfs_start_time);
  
  std::cout << "Finished BFS (" << bfs_time <<" seconds)"<< std::endl;

  validate_level(level, bfs_level_reference_file_name);
  std::cout << "Passed validation" << std::endl;

  utility::unmap_file(false,
		      utility::get_file_size(graph_file_name),
		      const_cast<uint64_t *>(index));
  return 0;
}
