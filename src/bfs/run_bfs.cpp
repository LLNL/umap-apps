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

#include "bfs_kernel.hpp"
#include "../utility/commandline.hpp"
#include "../utility/bitmap.hpp"
#include "../utility/umap_file.hpp"
#include "../utility/time.hpp"

void bfs_usage(char* pname)
{
  std::cerr << "BFS usage: " << pname << " -n <num vertices> -m <num edges>" << std::endl;
}

struct bfs_opts {
  size_t num_vertices;
  size_t num_edges;
};

int parse_options(void* optstruct, int argc, char **argv)
{
  bfs_opts* opts = (bfs_opts*) optstruct;

  int c;
  while ((c = getopt(argc, argv, "n:m:g:h")) != -1) {
    switch (c) {
      case 'n': /// Required
        opts->num_vertices = std::stoull(optarg);
        break;

      case 'm': /// Required
        opts->num_edges = std::stoull(optarg);
        break;

      case 'h':
        return -1;
        break;

      default:
        return -1;
    }
  }

  return 0;
}

std::pair<uint64_t *, uint64_t *>
map_graph(const size_t num_vertices, const size_t num_edges, const std::string &graph_file_name, bool usemmap) {
  const size_t graph_size = (num_vertices + 1 + num_edges) * sizeof(uint64_t);

  int fd = -1;
  void *map_raw_address = nullptr;
  map_raw_address = utility::map_in_file(graph_file_name, false, true, usemmap, graph_size);
  if (map_raw_address == nullptr) {
    std::cerr << "Failed to map the graph" << std::endl;
    std::abort();
  }

  uint64_t *index = static_cast<uint64_t *>(map_raw_address);
  const std::ptrdiff_t edges_offset = num_vertices + 1;
  uint64_t *edges = static_cast<uint64_t *>(map_raw_address) + edges_offset;

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
  utility::umt_optstruct_t global_opts;
  bfs_opts test_opts;
  utility::umt_getoptions(&global_opts, argc, argv);
  utility::umt_handle_options(&test_opts, &parse_options, &bfs_usage);

  const uint64_t *index = nullptr;
  const uint64_t *edges = nullptr;
  std::tie(index, edges) = map_graph(test_opts.num_vertices, test_opts.num_edges, std::string(global_opts.filename), global_opts.usemmap);

  std::vector<uint16_t> level(test_opts.num_vertices); // Array to store each vertex's level (a distance from the source vertex)
  std::vector<uint64_t> visited_filter(utility::bitmap_size(test_opts.num_vertices)); // bitmap data to store 'visited' information

  bfs::init_bfs(test_opts.num_vertices, level.data(), visited_filter.data());
  find_bfs_root(test_opts.num_vertices, index, level.data());

  const auto bfs_start_time = utility::elapsed_time_sec();
  const uint16_t max_level = bfs::run_bfs(test_opts.num_vertices, index, edges, level.data(), visited_filter.data());
  const auto bfs_time = utility::elapsed_time_sec(bfs_start_time);
  std::cout << "BFS took (s) " << bfs_time << std::endl;

  count_level(test_opts.num_vertices, max_level, level.data());

  return 0;
}
