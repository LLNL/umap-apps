/* This file is part of UMAP.  For copyright information see the COPYRIGHT
 * file in the top level directory, or at
 * https://github.com/LLNL/umap/blob/master/COPYRIGHT This program is free
 * software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.  This program is distributed in
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the terms and conditions of the GNU Lesser General Public License for
 * more details.  You should have received a copy of the GNU Lesser General
 * Public License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <atomic>
#include <fstream>
#include <cstring>

#include "umap/umap.h"
#include "umap/store/SparseStore.h"
#include "umap_file.hpp"

namespace utility {
  namespace{
    Umap::SparseStore* store;
  }

/* Remove the directory with all files under it */
int remove_directory(std::string root_path){
  DIR *directory;
  struct dirent *ent;
  if ((directory = opendir(root_path.c_str())) != NULL){
    while ((ent = readdir(directory)) != NULL){
      std::string filename_str(ent->d_name);
      std::string file_path = root_path + "/" + filename_str; 
      unlink(file_path.c_str());
    }
    closedir(directory);
    if( remove(root_path.c_str()) != 0 ){
      perror("Error deleting directory from new function");
      return -1;
    }
    else{
      std::cout << "Successfully deleted directory" << std::endl;
      return 0;
    }
  } else {
    perror("Error deleting directory");
    return -1; 
  }
}

/*
 * Wrapper function to map a region using umap with SparseStore
 */
  
  void* map_in_sparse_store(
      std::string root_path,
      bool initonly,
      bool noinit,
      uint64_t numbytes,
      void* start_addr,
      size_t file_size){

     void * region = NULL;
     
     if ( initonly || !noinit ) {
       struct stat info;
       if(stat( root_path.c_str(), &info ) == 0){
         remove_directory(root_path.c_str());
       }  
     }
    
     // Umap::SparseStore* store;
     size_t page_size = utility::get_umap_page_size();
     store = new Umap::SparseStore(numbytes,page_size,root_path,file_size);

     // Check status to make sure that the store object was able to open the directory
     if (store->get_directory_creation_status() != 0){
       std::cerr << "Error: Failed to create directory at " << root_path << std::endl;
       return NULL;
     }

     // call umap
     int flags = UMAP_PRIVATE;
     
     if (start_addr != nullptr)
      flags |= MAP_FIXED;

     const int prot = PROT_READ|PROT_WRITE;
     /* Here, the file descriptor handed top umap is -1, as we do not start with mapping a file
        instead, file(s) will be created incrementally as needed using the "store" object. */

     region = umap_ex(start_addr, numbytes, prot, flags, -1, 0, store);
     if ( region == UMAP_FAILED ) {
       std::ostringstream ss;
       ss << "umap_mf of " << numbytes
          << " bytes failed for " << root_path << ": ";
       perror(ss.str().c_str());
       delete store;
       return NULL;
     }

     return region;
   }

   void unmap_sparse_store(uint64_t numbytes, void* region){
     if (uunmap(region, numbytes) < 0) {
       std::ostringstream ss;
       ss << "uunmap of failure: ";
       perror(ss.str().c_str());
       exit(-1);
     }
     int sparse_store_close_files = store->close_files();
     if (sparse_store_close_files != 0 ){
       std::cerr << "Error closing SparseStore files" << std::endl;
       delete store;
       exit(-1);
     }
     delete store;
   }



}
