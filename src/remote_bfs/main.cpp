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
#include <cstring>
#include <cassert>
#include "mpi.h"

#include "util.hpp"
#include "umap/umap.h"
#include "umap/store/StoreNetwork.h"

#include "../utility/umap_file.hpp"
#include "../utility/time.hpp"
#include "../utility/file.hpp"

Umap::Store* ds;
int rank, num_proc;

size_t get_aligned_size( size_t original_size ){
  
  /* round up to page aligned size */
  size_t umap_page_size = umapcfg_get_umap_page_size();
  size_t aligned_pages  = (original_size - 1)/umap_page_size + 1;
  size_t aligned_size   = umap_page_size * aligned_pages;

  return aligned_size;
}


void create_datastore_server( std::string filename)
{

  int fd;
  int o_opts = O_RDONLY | O_LARGEFILE | O_DIRECT;
  if ( ( fd = open(filename.c_str(), o_opts, S_IRUSR | S_IWUSR) ) == -1 ) {
    std::string estr = "Failed to open " + filename + ": ";
    perror(estr.c_str());
  }

  struct stat sbuf;
  if (fstat(fd, &sbuf) == -1) {
    std::string estr = "Failed to get status (fstat) for " + filename + ": ";
    perror(estr.c_str());
  }

  off_t numbytes = (off_t)sbuf.st_size;

  const int prot = PROT_READ;
  void* region = NULL;
  int flags = MAP_SHARED | MAP_NORESERVE;
  region = mmap(NULL, numbytes, prot, flags, fd, 0);
  if (region == MAP_FAILED) {
    std::ostringstream ss;
    ss << "mmap of " << numbytes << " bytes failed for " << filename << ": ";
    perror(ss.str().c_str());
  }
  

  /* Create the network-based datastore */

  size_t umap_page_size = umapcfg_get_umap_page_size();
  size_t total_aligned_pages  = (numbytes - 1)/umap_page_size + 1;
  size_t pages_per_server = total_aligned_pages/num_proc;

  size_t page_st = pages_per_server*rank;
  size_t page_end = page_st + pages_per_server;
  if( rank==(num_proc-1) )
    page_end = total_aligned_pages;

  size_t aligned_size = (page_end-page_st+1)*umap_page_size;
  void* region_dup = malloc(aligned_size);
  if( !region_dup ){
    std::string estr = "Failed to allocate region_dup ";
    perror(estr.c_str());
  }

  std::cout << "Graph file size  " << (numbytes) <<" aligned to "
	    << total_aligned_pages<<" pages, "
	    << (total_aligned_pages*umap_page_size) << std::endl;

  /* Each Server fetch in its portion from the file system */
  size_t copy_bytes = aligned_size;
  if( rank==(num_proc-1) )
    copy_bytes = numbytes - aligned_size*rank;
  
  if( copy_bytes>0 ){
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int num_threads = omp_get_num_threads();
      size_t thread_stride = copy_bytes/num_threads;
      size_t copy_length = thread_stride;
      if( tid==(num_threads-1) )
	copy_length = copy_bytes - thread_stride*tid;
            
      memcpy( (char*)region_dup+thread_stride*tid,
	      (char*)region + rank*pages_per_server*umap_page_size + thread_stride*tid,
	      copy_length);

      if(tid==0)
	std::cout << num_threads << " threads finished memcpy" << std::endl;
    }
  }//End of copy from file to memory

  ds  = new Umap::StoreNetworkServer("graph_ptr", region_dup, aligned_size);
  std::cout << "Server " <<rank << " graph_ptr is Registed " << std::endl;
}


std::pair<uint64_t *, uint64_t *> map_graph_client(const bfs_options &options)
{

  /* umap page aligned size */
  size_t original_size  = (options.num_vertices + 1 + options.num_edges) * sizeof(uint64_t);
  size_t aligned_size   = get_aligned_size(original_size);

  /* Create the network-based datastore */
  ds  = new Umap::StoreNetworkClient("graph_ptr", 0);
  size_t numbytes = ds->get_size();
  std::cout << " Client graph_ptr is Registed "<<numbytes <<" bytes (aligned to "<<aligned_size<<" )\n\n";
  assert(aligned_size<=numbytes);

  
  /* create a UMap region on the datastore */
  const int prot = PROT_READ;
  int flags = UMAP_PRIVATE;
  void *const map_raw_address = Umap::umap_ex(NULL, aligned_size, prot, flags, -1, 0, ds);
  if ( map_raw_address == UMAP_FAILED ) {
    std::ostringstream ss;
    ss << "umap " << options.graph_file_name << " from network-based datastore ";
    perror(ss.str().c_str());
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

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  
  if(options.is_server){

    create_datastore_server(options.graph_file_name);

    while(true)
      sleep(10);
    
  }
  else{

    disp_umap_env_variables();

    const uint64_t *index = nullptr;
    const uint64_t *edges = nullptr;
    std::tie(index, edges) = map_graph_client(options);

    // Array to store each vertex's level (a distance from the source vertex)
    std::vector<uint16_t> level(options.num_vertices);

    // bitmap data to store 'visited' information
    std::vector<uint64_t> visited_filter(utility::bitmap_size(options.num_vertices));
    
    init_bfs(options.num_vertices, level.data(), visited_filter.data());
    find_bfs_root(options.num_vertices, index, level.data());

    std::cout << "Before BFS #of page faults" << std::endl;
    print_num_page_faults();
    const auto bfs_start_time = utility::elapsed_time_sec();
    const uint16_t max_level = run_bfs(options.num_vertices,
				       index,
				       edges,
				       level.data(),
				       visited_filter.data());
    const auto bfs_time = utility::elapsed_time_sec(bfs_start_time);
    std::cout << "Client "<< rank <<" BFS took (s)\t" << bfs_time << std::endl;
    std::cout << "After BFS #of page faults" << std::endl;
    print_num_page_faults();
  
    /*Start validation */
    count_level(options.num_vertices, max_level, level.data());
    MPI_Barrier(MPI_COMM_WORLD);

    /* Unmap file */
    if ( uunmap( (void*)index, 0) ) {
      int eno = errno;
      std::cerr << "Failed to unmap network datastore: " << strerror(eno) << std::endl;
      return -1;
    }

  }

  delete(ds);

  return 0;
  
}


