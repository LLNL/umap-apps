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
#include <list>
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
#include <pthread.h>
#include "umap/umap.h"
#include "fitsio.h"

#include "../utility/commandline.hpp"
#include "umap/Store.h"

namespace utility {
namespace umap_fits_file {

class CfitsStoreFile;

struct Point {
  Point(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}
  Point(int64_t x, int64_t y) : x(x), y(y), z(0) {}
  Point(const Point& p) : x(p.x), y(p.y), z(p.z) {}
  Point() : x(0), y(0), z(0) {}
  int64_t x, y, z;
  bool operator>(const Point& lhs) const {
    return (x > lhs.x || y > lhs.y || z > lhs.z);
  }
  Point& operator=(const Point& lhs) {
    this->x = lhs.x; this->y = lhs.y; this->z = lhs.z;
    return *this;
  }
};

struct Fits_Base {
  Point dim; // the pixel dimensions of the object's memory space
  Point max; // the pixel dimensions of the bounding box for the image data in the object
  Point pos; // the absolute position of the object (in pixels)
  std::size_t elem_size; // the size in bytes of each pixel
};

struct Tile_File {
  int fd;
  std::string fname;
  std::size_t tile_start;
  std::size_t tile_size;
};

class Tile : public Fits_Base {
friend std::ostream &operator<<(std::ostream &os, utility::umap_fits_file::Tile const &ft);
friend class Layer;
public:
  Tile(const std::string& _fn);
  ssize_t buffered_read(void*, std::size_t, off_t);
  ssize_t read_compressed(void*, std::size_t, off_t);
private:
  Tile_File file;
  void* map;
  std::size_t map_start; // start of the file data in the map
  std::size_t map_size; // size of the map
  std::string rawdate;
  double psf;
  double exptime;
  bool compressed;
  long comp_readlen;
  fitsfile* fptr;
  pthread_mutex_t fits_mtx;
};
std::ostream &operator<<(std::ostream &os, utility::umap_fits_file::Tile const &ft);

class Layer : public Fits_Base {
  friend class CfitsStoreFile;
public:
  Layer() { this->page_size = utility::umt_getpagesize(); }
  void add_tile(const std::string& _fn);
  ssize_t buffered_read(void* request_buf, std::size_t request_size, off_t request_offset);
private:
  std::list<utility::umap_fits_file::Tile> tiles;
  off_t page_size;
};

struct Cube : public Fits_Base {
  Cube() : layer_size(0), cube_size(0) {}
  size_t layer_size;  // Size of each layer (assumed to be the same for each layer)
  size_t cube_size;  // Total bytes in cube
  off_t page_size;
  vector<utility::umap_fits_file::Layer> layers;  // Just one column for now
};

static std::unordered_map<void*, Cube*>  Cubes;

class CfitsStoreFile : public Store {
  public:
    CfitsStoreFile(Cube* _cube_, size_t _rsize_, size_t _aligned_size)
      : cube{_cube_}, rsize{_rsize_}, aligned_size{_aligned_size}{}

    ssize_t read_from_store(char* buf, size_t nb, off_t off) {
      ssize_t rval = 0;
      off_t layerno = off / cube->layer_size;
      off_t layeroffset = off % cube->layer_size;
      
//       std::cout << off << ":\t" << " t:" << layerno 
//       << "\tx:" << ((layeroffset / cube->layers[layerno].elem_size) % (cube->layers[layerno].dim.x))
//       << "\ty:" << ((layeroffset / cube->layers[layerno].elem_size) / (cube->layers[layerno].dim.x))
//       << std::endl;

      if ( ( rval = cube->layers[layerno].buffered_read(buf, nb, layeroffset) ) == -1) {
        perror("ERROR: buffered_read failed");
        exit(1);
      }

      // Fill the rest with NaNs if the read returned less than the requested amount.
      if (rval < nb) {
        memset((void*)&buf[rval], 0xff, nb - rval - 1);
      }
      return rval;
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

  Cube* cube = new Cube();
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
    
    cube->layers.emplace(cube->layers.end());
    Layer& L = cube->layers.back();
    L.add_tile(ss.str());
    ss.str(""); ss.clear();

    utility::umap_fits_file::Point dim = L.dim;
    if ( *BytesPerElement == 0 ) {
      *xDim = dim.x;
      *yDim = dim.y;
      *BytesPerElement = L.elem_size;
      cube->layer_size = (dim.x * dim.y * L.elem_size);
      cube->max = L.max;
      cube->dim = L.dim;
    }
    else {
      assert( *xDim == dim.x && *yDim == dim.y && *BytesPerElement == L.elem_size );
      assert( (dim.x * dim.y * L.elem_size) == cube->layer_size );
      if (L.max > cube->max) {
        cube->max = L.max;
        cube->dim = L.dim;
      }
    }
    *zDim = i;

    cube->cube_size += cube->layer_size;
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

void Layer::add_tile(const std::string& _fn) {
  this->tiles.emplace(this->tiles.end(), _fn);
  const Tile& T = this->tiles.back();
  if (T.max > this->max) {
    this->dim = T.dim;
    this->max = T.max;
    this->elem_size = T.elem_size;
    // align the layer size to the umap page size to reduce address ambiguity
    off_t pgaligned = ((this->dim.x * this->elem_size) & ~(this->page_size-1));
    this->dim.x = (pgaligned + this->page_size) / this->elem_size;
  }
}

Tile::Tile(const std::string& _fn) : fits_mtx(PTHREAD_MUTEX_INITIALIZER)
{
  int status = 0;
  LONGLONG headstart;
  LONGLONG datastart;
  LONGLONG dataend;
  int bitpix;
  long naxis[2];
  int naxes;
  int hdutype;
  char buf[80];
  int open_flags = (O_RDONLY | O_LARGEFILE | O_DIRECT);

  file.fname = _fn;
  file.tile_start = (size_t)0;
  file.tile_size = (size_t)0;
  dim.x = (size_t)0;
  dim.y = (size_t)0;
  elem_size = 0;
  file.fd = -1;

  if ( fits_open_data(&fptr, file.fname.c_str(), READONLY, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_get_hduaddrll(fptr, &headstart, &datastart, &dataend, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if (dataend - datastart < 2880*16) {
    // data is too small; switch to extension 2 and continue reading headers
    if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) {
      fits_report_error(stderr, status);
      exit(-1);
    }
  }

  if ( fits_get_img_type(fptr, &bitpix, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_get_img_param(fptr, 2, &bitpix, &naxes, naxis, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_read_key_str(fptr, "DATE-AVG", buf, 0, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  } else {
    this->rawdate = std::string(buf);
  }
  if ( fits_read_key_dbl(fptr, "PSF_FWHM", &this->psf, 0, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  if ( fits_read_key_dbl(fptr, "EXPTIME", &this->exptime, 0, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_read_key_lnglng(fptr, "CRVAL1A", (long long int*) &this->pos.x, 0, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_read_key_lnglng(fptr, "CRVAL2A", (long long int*) &this->pos.y, 0, &status) ) {
    fits_report_error(stderr, status);
    exit(-1);
  }
  
  if ( !fits_read_key_str(fptr, "ZCMPTYPE", buf, 0, &status) ) {
    std::cerr << "Warning: " << file.fname << " is compressed." << std::endl;
    // file is compressed
    this->compressed = true;
    if ( fits_read_key_lnglng(fptr, "ZTILE1", (long long int*) &this->comp_readlen, 0, &status) ) {
      fits_report_error(stderr, status);
      exit(-1);
    }
  } else {
    this->compressed = false;
  }

  if ( ( file.fd = open(file.fname.c_str(), open_flags) ) == -1 ) {
    perror(file.fname.c_str());
    exit(-1);
  }

  dim.x = (size_t)naxis[0];
  dim.y = (size_t)naxis[1];
  elem_size = bitpix < 0 ? (size_t)( ( bitpix * -1 ) / 8 ) : (size_t)( bitpix / 8 );
  file.tile_start = (size_t)datastart;
  file.tile_size = (size_t)(dim.x * dim.y * elem_size);
  this->max = this->dim;

  std::size_t pgaligned_tile_start = file.tile_start & ~(4096-1);
  map_start = file.tile_start - pgaligned_tile_start;
  map_size = file.tile_size + map_start;
  
  if (!this->compressed) {
    if ( ( map = mmap(0, map_size, PROT_READ, MAP_PRIVATE | MAP_NORESERVE, file.fd, pgaligned_tile_start) ) == 0 ) {
      perror(file.fname.c_str());
      exit(-1);
    }
    // Set some advice flags to minimize unwanted buffering and read-ahead on
    // sparse data accesses
    posix_fadvise(file.fd, pgaligned_tile_start, map_size, POSIX_FADV_RANDOM);
    madvise(map, map_size, MADV_RANDOM | MADV_DONTDUMP);
    assert( (dataend - datastart) >= (dim.x * dim.y * elem_size) );
  }
}

ssize_t Layer::buffered_read(void* request_buf, std::size_t request_size, off_t request_offset)
{
  Tile& t = this->tiles.back(); // currently only one tile
  uint64_t t_row_size = t.dim.x * t.elem_size;
  uint64_t elem_offset = request_offset / t.elem_size;
  uint64_t req_row = (elem_offset / this->dim.x);
  uint64_t req_col = (elem_offset % this->dim.x);
  uint64_t req_size = request_size;
  assert(req_col < t.dim.x); // layer is page-aligned in x, so this should not happen.
  if ((req_col + request_size) > t_row_size) {
    req_size = t_row_size - req_col;
  }
  uint64_t t_req_offset = (req_row * t.dim.x + req_col) * t.elem_size;
  ssize_t rv = t.buffered_read(request_buf, req_size, t_req_offset);
  if (rv < request_size) {
    memset(&((char*)request_buf)[rv], 0xff, request_size - rv);
  }
}

ssize_t Tile::buffered_read(void* request_buf, std::size_t request_size, off_t request_offset)
{
  if (this->compressed) {
    return this->read_compressed(request_buf, request_size / elem_size, request_offset / elem_size);
  }
  std::size_t ua_request_size = request_size + map_start;
  off_t eof = 0;
  off_t req_end = request_offset + ua_request_size;

  // Make sure that the memcpy doesn't go past EOF
  if (req_end >= map_size) {
    eof = req_end - map_size;
  }

  void* request = &((char*)map)[request_offset];
  memcpy(request_buf, request, ua_request_size - eof);
  madvise(request, request_size*2, MADV_DONTNEED);
  // Realign the data in the request buffer.
  // TODO: Find a better way to fix this.
  // This is a hack that relies on the fact that copy_buf passed in by
  // the userfault handler from umap is actually allocated to be 2 umap pages
  // rather than the request_size (1 page) given to this function.
  memmove(request_buf, &((char*)request_buf)[map_start], request_size);

  return request_size - eof;
}

ssize_t Tile::read_compressed(void* request_buf, std::size_t request_size, off_t request_offset) {
  // Warning: compressed reads are not thread-safe
  int anynul = 0;
  int status = 0;
  long fpixel[2] = {
    (request_offset % this->dim.x) + 1,
    (request_offset / this->dim.x) + 1
  };
  long lpixel[2] = { fpixel[0] + (long)request_size, fpixel[1] };
  if (lpixel[0] > this->dim.x) {
    lpixel[0] = this->dim.x;
  }
  assert(lpixel[0] - fpixel[0] <= request_size);
  long rqsize = lpixel[0] - fpixel[0];
  pthread_mutex_lock( &fits_mtx );
  fits_read_pix(fptr, TFLOAT, fpixel, rqsize, 0, (float*)request_buf, &anynul, &status);
  pthread_mutex_unlock( &fits_mtx );
  return rqsize * this->elem_size;
}

std::ostream &operator<<(std::ostream &os, Tile const &ft)
{
  os << ft.file.fname << " "
     << "Start=" << ft.file.tile_start << ", "
     << "Size=" << ft.file.tile_size << ", "
     << "XDim=" << ft.dim.x << ", "
     << "YDim=" << ft.dim.y << ", "
     << "ESize=" << ft.elem_size << " ";

  return os;
}

}
}
#endif
