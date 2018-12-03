/*
 * This file is part of UMAP.
 *
 * For copyright information see the COPYRIGHT file in the top level
 * directory or at https://github.com/LLNL/umap/blob/master/COPYRIGHT.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License (as published by
 * the Free Software Foundation) version 2.1 dated February 1999.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307 USA
 */
#ifndef _UMAP_PERFITS_H
#define _UMAP_PERFITS_H
#include <ostream>
#include <string>
#include <vector>
#include <stdint.h>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include "umap/umap.h"
#include "fitsio.h"

#include "../utility/commandline.hpp"
#include "umap/Store.h"

namespace utility {
namespace umap_fits_file {

struct Tile_Dim {
  std::size_t xDim;
  std::size_t yDim;
  std::size_t elem_size;
};

struct Tile_File {
  int fd;
  std::string fname;
  std::size_t tile_start;
  std::size_t tile_size;
};

class Tile {
friend std::ostream &operator<<(std::ostream &os, utility::umap_fits_file::Tile const &ft);
public:
  Tile(const std::string& _fn);
  ssize_t buffered_read(char*, std::size_t, void*, std::size_t, off_t);
  Tile_Dim get_Dim() { return dim; }
private:
  Tile_File file;
  Tile_Dim  dim;
};
std::ostream &operator<<(std::ostream &os, utility::umap_fits_file::Tile const &ft);

struct Cube {
  size_t tile_size;  // Size of each tile (assumed to be the same for each tile)
  size_t cube_size;  // Total bytes in cube
  off_t page_size;
  vector<utility::umap_fits_file::Tile> tiles;  // Just one column for now
};

static std::unordered_map<void*, Cube*>  Cubes;

class CfitsStoreFile : public Store {
  public:
    CfitsStoreFile(Cube* _cube_, size_t _rsize_, size_t _aligned_size)
      : cube{_cube_}, rsize{_rsize_}, aligned_size{_aligned_size}
    {
      if ( posix_memalign((void**)(&aligned_buf), aligned_size, aligned_size) ) {
            exit(1);
      }
    }

    ~CfitsStoreFile() {
      free(aligned_buf);
    }

    ssize_t read_from_store(char* buf, size_t nb, off_t off) {
      ssize_t rval;
      ssize_t bytesread = 0;

      //
      // Now read in remaining bytes
      //
      while ( nb ) {
        off_t tileno = off / cube->tile_size;
        off_t tileoffset = off % cube->tile_size;
        size_t bytes_to_eof = cube->tile_size - tileoffset;
        size_t bytes_to_read = std::min(bytes_to_eof, nb);

        if ( ( rval = cube->tiles[tileno].buffered_read(aligned_buf, aligned_size, buf, bytes_to_read, tileoffset) ) == -1) {
          perror("ERROR: buffered_read failed");
          exit(1);
        }

        bytesread += rval;
        nb -= rval;
        buf += rval;
        off += rval;
      }

      return bytesread;
    }

    ssize_t  write_to_store(char* buf, size_t nb, off_t off) {
      assert("FITS write not supported" && 0);
      return 0;
    }

    void* region;

  private:
    Cube* cube;
    char* aligned_buf;
    size_t rsize;
    size_t aligned_size;
    int fd;
};

/* Returns pointer to cube[Z][Y][X] Z=time, X/Y=2D space coordinates */
void* PerFits_alloc_cube(
    string name,
    size_t* BytesPerElement,            /* Output: size of each element of cube */
    size_t* xDim,                       /* Output: Dimension of X */
    size_t* yDim,                       /* Output: Dimension of Y */
    size_t* zDim                        /* Output: Dimension of Z */
)
{
  void* region = NULL;

  Cube* cube = new Cube{.tile_size = 0, .cube_size = 0};
  cube->page_size = utility::umt_getpagesize();
  string basename(name);

  *xDim = *yDim = *BytesPerElement = 0;
  for (int i = 1; ; ++i) {
    std::stringstream ss;
    ss << basename << i << ".fits";
    struct stat sbuf;

    if ( stat(ss.str().c_str(), &sbuf) == -1 ) {
      if ( i == 1 ) {
        cerr << "File: " << ss.str() << " does not exist\n";
        return region;
      }
      break;
    }
    utility::umap_fits_file::Tile T( ss.str() );
    ss.str(""); ss.clear();

    ss << T << endl;

    utility::umap_fits_file::Tile_Dim dim = T.get_Dim();
    if ( *BytesPerElement == 0 ) {
      *xDim = dim.xDim;
      *yDim = dim.yDim;
      *BytesPerElement = dim.elem_size;
      cube->tile_size = (dim.xDim * dim.yDim * dim.elem_size);
    }
    else {
      assert( *xDim == dim.xDim && *yDim == dim.yDim && *BytesPerElement == dim.elem_size );
    }
    *zDim = i;

    cube->cube_size += cube->tile_size;

    cube->tiles.push_back(T);
  }

  // Make sure that our cube is padded if necessary to be page aligned
  
  size_t psize = utility::umt_getpagesize();
  long remainder = cube->cube_size % psize;

  cube->cube_size += remainder ? (psize - remainder) : 0;


  CfitsStoreFile* cstore;
  cstore = new CfitsStoreFile{cube, cube->cube_size, psize};

  const int prot = PROT_READ|PROT_WRITE;
  int flags = UMAP_PRIVATE;

  cstore->region = umap_ex(NULL, cube->cube_size, prot, flags, 0, 0, cstore);
  if ( cstore->region == UMAP_FAILED ) {
      ostringstream ss;
      ss << "umap of " << cube->cube_size << " bytes failed for Cube";
      perror(ss.str().c_str());
      return NULL;
  }

  Cubes[cstore->region] = cube;
  return cstore->region;
}

void PerFits_free_cube(void* region)
{
  auto it = Cubes.find(region);
  assert( "free_cube: failed to find control object" && it != Cubes.end() );
  Cube* cube = it->second;

  if (uunmap(region, cube->cube_size) < 0) {
    ostringstream ss;
    ss << "uunmap of " << cube->cube_size << " bytes failed on region " << region << ": ";
    perror(ss.str().c_str());
    exit(-1);
  }

  delete cube;

  Cubes.erase(region);
}

Tile::Tile(const std::string& _fn)
{
  fitsfile* fptr = NULL;
  int status = 0;
  LONGLONG headstart;
  LONGLONG datastart;
  LONGLONG dataend;
  int bitpix;
  long naxis[2];
  int naxes;
  int open_flags = (O_RDONLY | O_LARGEFILE | O_DIRECT);

  file.fname = _fn;
  file.tile_start = (size_t)0;
  file.tile_size = (size_t)0;
  dim.xDim = (size_t)0;
  dim.yDim = (size_t)0;
  dim.elem_size = 0;
  file.fd = -1;

  if ( fits_open_data(&fptr, file.fname.c_str(), READONLY, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_get_hduaddrll(fptr, &headstart, &datastart, &dataend, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_get_img_type(fptr, &bitpix, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_get_img_param(fptr, 2, &bitpix, &naxes, &naxis[0], &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_close_file(fptr, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( ( file.fd = open(file.fname.c_str(), open_flags) ) == -1 ) {
    perror(file.fname.c_str());
    exit(-1);
  }

  dim.xDim = (size_t)naxis[0];
  dim.yDim = (size_t)naxis[1];
  dim.elem_size = bitpix < 0 ? (size_t)( ( bitpix * -1 ) / 8 ) : (size_t)( bitpix / 8 );
  file.tile_start = (size_t)datastart;
  file.tile_size = (size_t)(dim.xDim * dim.yDim * dim.elem_size);

  assert( (dataend - datastart) >= (dim.xDim * dim.yDim * dim.elem_size) );
}

ssize_t Tile::buffered_read(char* aligned_buf, std::size_t aligned_size, void* request_buf, std::size_t request_size, off_t request_offset)
{
  off_t data_offset = request_offset + file.tile_start;
  off_t aligned_data_offset = data_offset & ~(aligned_size - 1);
  off_t copy_start_offset = data_offset - aligned_data_offset;
  off_t copy_amount = aligned_size - copy_start_offset;

  assert("buffered_read: invalid offset received" 
      && aligned_data_offset >= 0 
      && aligned_data_offset < (file.tile_start + file.tile_size));

  off_t bytes_read;
  if ((bytes_read = ::pread(file.fd, aligned_buf, aligned_size, aligned_data_offset)) == -1) {
    perror("ERROR: pread failed");
    exit(1);
  }

  assert("Incorrect number of bytes read" && copy_amount > 0);

  copy_amount = std::min(bytes_read, copy_amount);
  copy_amount = std::min((off_t)request_size, copy_amount);

  memcpy(request_buf, &aligned_buf[copy_start_offset], copy_amount);

  return copy_amount;
}

std::ostream &operator<<(std::ostream &os, Tile const &ft)
{
  os << ft.file.fname << " "
     << "Start=" << ft.file.tile_start << ", "
     << "Size=" << ft.file.tile_size << ", "
     << "XDim=" << ft.dim.xDim << ", "
     << "YDim=" << ft.dim.yDim << ", "
     << "ESize=" << ft.dim.elem_size << " ";

  return os;
}
}
}
#endif
