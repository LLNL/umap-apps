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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <cassert>
#include <utility>

# include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

#include "umap/umap.h"
#include "fitsio.h"

#include "../utility/commandline.hpp"
#include "umap/Store.h"

namespace utility { namespace umap_fits_file_internal {
// internal interface to umap_fits_file - contains umap handler
class CfitsStoreFile;
struct Cube;

// Timestamp with an epoch at 2000-01-01 00:00:00.000
struct fits_time {
private:
  int64_t t_abs;
  uint64_t Y, M, D, h, m, s, ns;
public:
  fits_time() = default;
  fits_time(const fits_time& c) = default;
  void parse(std::string t) {
    std::sscanf(t.c_str(), "%lu-%lu-%luT%lu:%lu:%lu.%lu", &Y, &M, &D, &h, &m, &s, &ns);
    this->t_abs = this->toint();
  }
  fits_time operator- (const fits_time& rhs) const {
    fits_time rop(*this);
    rop.Y = this->Y - rhs.Y;
    rop.M = this->M - rhs.M;
    rop.D = this->D - rhs.D;
    rop.h = this->h - rhs.h;
    rop.m = this->m - rhs.m;
    rop.s = this->s - rhs.s;
    rop.ns = this->ns - rhs.ns;
    rop.t_abs = rop.toint();
    return rop;
  }
  int64_t toint() const {
    int64_t months = (this->Y - 2000) * 12 + this->M;
    int64_t days   = months * 30 + this->D;
    int64_t hours  = days * 24 + this->h;
    int64_t mins   = hours * 60 + this->m;
    int64_t secs   = mins * 60 + this->s;
    int64_t nsecs  = secs * 1000000000 + this->ns;
    return nsecs / 1000000; // msecs
  }
  
  bool operator<(const fits_time& rhs) const {
    return this->t_abs < rhs.t_abs;
  }
  bool operator==(const fits_time& rhs) const {
    return this->t_abs == rhs.t_abs;
  }
};

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
  Point operator+(const Point& lhs) const {
    Point r(*this);
    r.x += lhs.x; r.y += lhs.y; r.z += lhs.z;
    return r;
  }
};

struct Fits_Base {
  Fits_Base() : elem_size(0) {}
  Point dim; // the pixel dimensions of the object's memory space
  Point max; // the pixel dimensions of the bounding box for the image data in the object
  Point pos; // the absolute position of the object (in pixels)
  std::size_t elem_size; // the size in bytes of each pixel
};

struct Layer_Base : public Fits_Base {
  Layer_Base() {};
  fits_time ts;
  bool operator<(const Layer_Base& rhs) const {
    return this->ts < rhs.ts;
  }
  bool operator==(const Layer_Base& rhs) const {
    return this->ts == rhs.ts;
  }
  
};

struct Tile_File {
  int fd;
  std::string fname;
  std::size_t tile_start;
  std::size_t tile_size;
};

class Tile : public Layer_Base {
friend std::ostream &operator<<(std::ostream &os, utility::umap_fits_file_internal::Tile const &ft);
friend class Layer;
friend struct Cube;
public:
  Tile() = default;
  Tile(const std::string& _fn);
  std::size_t read_row(char*, std::size_t, std::size_t, std::size_t) const;
  std::size_t read_row_compressed(char*, std::size_t, std::size_t, std::size_t) const;
private:
  Tile_File file;
  void* map;
  std::size_t data_offset; // offset of the file data
  std::size_t map_start;   // start of the file data in the map
  std::size_t map_size;    // size of the map
  double psf;
  double exptime;
  bool compressed;
  long comp_readlen;
  fitsfile* fptr;
  pthread_mutex_t fits_mtx;
};
std::ostream &operator<<(std::ostream &os, utility::umap_fits_file_internal::Tile const &ft);

class Layer : public Layer_Base {
  friend class CfitsStoreFile;
  friend struct Cube;
public:
  Layer(fits_time ts) {
    this->page_size = utility::umt_getpagesize();
    this->ts = ts;
  }

  void refit();
  void add_tile(Tile& t);
  const Tile* resolve_coords(std::size_t lx, std::size_t ly) const;
  const Tile* get_rnd_tile(std::mt19937 rnd_engine) const;
  std::size_t buffered_read(void* request_buf, std::size_t request_size, off_t request_offset) const;
private:
  std::vector<utility::umap_fits_file_internal::Tile*> tiles;
  off_t page_size;
};

}} // utility::umap_fits_file_internal

namespace std {
template<> struct hash<utility::umap_fits_file_internal::Layer> {
    size_t operator()(utility::umap_fits_file_internal::Layer_Base const& s) const noexcept {
      return s.ts.toint();
    }
};
template<> struct less<utility::umap_fits_file_internal::Layer> {
  using value_type = utility::umap_fits_file_internal::Layer_Base;
    size_t operator()( value_type const& s1, value_type const& s2) const noexcept {
      return s1 < s2;
    }
};
template<> struct hash<utility::umap_fits_file_internal::Tile> {
    size_t operator()(utility::umap_fits_file_internal::Layer_Base const& s) const noexcept {
      return s.ts.toint();
    }
};
template<> struct less<utility::umap_fits_file_internal::Tile> {
  using value_type = utility::umap_fits_file_internal::Layer_Base;
    size_t operator()( value_type const& s1, value_type const& s2) const noexcept {
      return s1 < s2;
    }
};
}

namespace utility { namespace umap_fits_file_internal {

struct Cube : public Fits_Base {
  friend class CfitsStoreFile;
public:
  Cube() : layer_size(0), cube_size(0) {}
  void add_tile(const std::string& _fn);
  void refit();
  const Tile* resolve_coords( std::size_t gx, std::size_t gy, std::size_t gz ) const;
  uint64_t index(Point xyz) const; // translate xyz to pixel offset
  std::tuple<std::size_t, std::size_t, std::size_t> get_rnd_coord(std::mt19937 rnd_engine) const;
  std::size_t layer_size;  // Size of each layer (in bytes, assumed to be the same for each layer)
  std::size_t cube_size;  // Total bytes in cube
  std::size_t page_size;
private:
  std::unordered_set<Layer> layer_set;
  std::list<Tile> tile_alloc;
  std::vector<const Layer*> layers;
  std::vector<Tile*> tiles;
};

class CfitsStoreFile : public Store {
  public:
    CfitsStoreFile(Cube* _cube_, std::size_t _rsize_, std::size_t _aligned_size)
      : cube{_cube_}, rsize{_rsize_}, aligned_size{_aligned_size}{}

    ssize_t read_from_store(char* buf, size_t nb, off_t off) {
      ssize_t rval = 0;
      off_t layerno = off / cube->layer_size;
      off_t layeroffset = off % cube->layer_size;
      
//       std::cout << off << ":\t" << " t:" << layerno 
//       << "\tx:" << ((layeroffset / cube->layers[layerno].elem_size) % (cube->layers[layerno].dim.x))
//       << "\ty:" << ((layeroffset / cube->layers[layerno].elem_size) / (cube->layers[layerno].dim.x))
//       << std::endl;

      if ( ( rval = cube->layers[layerno]->buffered_read(buf, nb, layeroffset) ) == -1) {
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

uint64_t Cube::index(Point xyz) const {
  uint64_t rv = (((xyz.z * this->dim.y) + xyz.y) * this->dim.x) + xyz.x;
  return rv;
}

void Cube::refit() {
  const int64_t nmax = std::numeric_limits<long long>::max();
  Point max(0, 0, this->layer_set.size()), min(nmax, nmax, 0), pmax(0, 0), pmin(nmax, nmax);
//   std::cout <<  this->layer_set.size() <<  std::endl;
  // layers only set pos and max when created.
  this->elem_size = this->layer_set.begin()->elem_size;
  for (auto it = this->tiles.begin(); it != this->tiles.end(); it++) {
    Tile* T = *it;
    if      ((T->pos.x + T->max.x) > max.x) { max.x = (T->pos.x + T->max.x); }
    else if ((T->pos.x + T->max.x) < min.x) { min.x = (T->pos.x + T->max.x); }
    if      ((T->pos.y + T->max.y) > max.y) { max.y = (T->pos.y + T->max.y); }
    else if ((T->pos.y + T->max.y) < min.y) { min.y = (T->pos.y + T->max.y); }
      
    if      (T->pos.x > pmax.x) { pmax.x = T->pos.x; }
    else if (T->pos.x < pmin.x) { pmin.x = T->pos.x; }
    if      (T->pos.y > pmax.y) { pmax.y = T->pos.y; }
    else if (T->pos.y < pmin.y) { pmin.y = T->pos.y; }
      
    assert(T->elem_size == this->elem_size);
  }
  
  this->pos.x = pmin.x; this->pos.y = pmin.y;
  this->max.x = max.x - this->pos.x; this->max.y = max.y - this->pos.y; 
  this->dim = this->max;
//   std::cout << min.x << " " << min.y << " " << min.z << std::endl;
//   std::cout << max.x << " " << max.y << " " << max.z << std::endl;
//   std::cout << this->max.x << " " << this->max.y << " " << this->max.z << " " << this->elem_size << std::endl;
//   std::cout << this->dim.x << " " << this->dim.y << " " << this->dim.z << " " << this->elem_size << std::endl;

  // pad the x dimension to a multiple of the umap page size
  this->max.x = this->dim.x;
  this->max.y = this->dim.y;
  this->max.z = this->dim.z = max.z;
  this->dim.x = ((this->dim.x * this->elem_size) & ~(this->page_size-1));
  if (this->dim.x & (this->page_size-1)) {
    this->dim.x += page_size;
  }
  this->dim.x /= this->elem_size;

  // make layer coordinates relative to cube
  for (auto it = this->layer_set.begin(); it != this->layer_set.end(); it++) {
    // this doesn't change the hash value. strict constness is unnecessary...
    Layer* l = ((Layer*)((void*) &*it));
    l->dim.x = dim.x; l->max.x = max.x; l->pos.x = 0;
    l->dim.y = dim.y; l->max.y = max.y; l->pos.y = 0;
    // make tile coordinates relative to layer coordinates
    for (auto it = l->tiles.begin(); it != l->tiles.end(); it++) {
      Tile* t = *it;
      t->pos.x -= l->pos.x; t->pos.y -= l->pos.y; t->pos.z -= l->pos.z;
    }
    this->layers.push_back(&(*it));
  }
  this->layer_size = this->dim.x * this->dim.y * this->elem_size;
  this->cube_size = this->layer_size * this->dim.z;
}

void Cube::add_tile(const std::string& _fn) {
  auto it = this->tile_alloc.emplace(this->tile_alloc.end(), _fn);
  Tile& T = *it;
  this->tiles.push_back(&T);
  auto r = this->layer_set.emplace(T.ts);
  auto dst_layer = r.first;

  // this doesn't change the hash value. cast away unnecessary unordered_set constness.
  Layer* dl = (Layer*)((void*) &(*dst_layer));
  dl->add_tile(T);
}

void Layer::add_tile(Tile& T) {
  this->tiles.emplace_back(&T);
  if (!this->elem_size) {
    this->elem_size = T.elem_size;
  } else {
    assert(this->elem_size == T.elem_size);
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
    std::cerr << "DATE-AVG" <<  std::endl;
    fits_report_error(stderr, status);
    exit(-1);
  } else {
    this->ts.parse(std::string(buf));
//     std::cout << buf << "\t" << this->ts.toint() << std::endl;
  }
  if ( fits_read_key_dbl(fptr, "PSF_FWHM", &this->psf, 0, &status) ) {
    std::cerr << "PSF_FWHM" <<  std::endl;
    fits_report_error(stderr, status);
    exit(-1);
  }
  if ( fits_read_key_dbl(fptr, "EXPTIME", &this->exptime, 0, &status) ) {
    std::cerr << "EXPTIME" <<  std::endl;
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_read_key_lnglng(fptr, "CRVAL1A", (long long int*) &this->pos.x, 0, &status) ) {
    std::cerr << "CRVAL1A" <<  std::endl;
    fits_report_error(stderr, status);
    exit(-1);
  }

  if ( fits_read_key_lnglng(fptr, "CRVAL2A", (long long int*) &this->pos.y, 0, &status) ) {
    std::cerr << "CRVAL2A" <<  std::endl;
    fits_report_error(stderr, status);
    exit(-1);
  }
  
  if ( !fits_read_key_str(fptr, "ZCMPTYPE", buf, 0, &status) ) {
//     std::cerr << "Warning: " << file.fname << " is compressed." << std::endl;
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
    close(file.fd);
    fits_close_file(fptr, &status);
    fptr = 0;
  }
}

const Tile* Cube::resolve_coords( std::size_t gx, std::size_t gy, std::size_t gz ) const {
  if ( gz < this->max.z ) {
    return this->layers[gz]->resolve_coords(gx, gy);
  } else {
    return 0;
  }
}

const Tile* Layer::resolve_coords(std::size_t gx, std::size_t gy) const {
  for (auto it = this->tiles.begin(); it != this->tiles.end(); it++) {
    const Tile& t = **it;
    if (
      gx >= t.pos.x && gy >= t.pos.y &&
      gx < (t.pos.x + t.max.x) && 
      gy < (t.pos.y + t.max.y)
    ) {
      return &t;
    }
  }
  return 0;
}

std::tuple<std::size_t, std::size_t, std::size_t> Cube::get_rnd_coord(std::mt19937 rnd_engine) const {
  std::uniform_int_distribution<std::size_t> layer_dist(0, this->layers.size() - 1);

  std::size_t pt_k  = layer_dist(rnd_engine);
  const Layer* L = this->layers[pt_k];
  int i = 0;

  const Tile* T = L->get_rnd_tile(rnd_engine);
  std::uniform_int_distribution<std::size_t> x_dist(T->pos.x, T->pos.x+T->max.x - 1);
  std::uniform_int_distribution<std::size_t> y_dist(T->pos.y, T->pos.y+T->max.y - 1);

  std::size_t pt_y = y_dist(rnd_engine);
  std::size_t pt_x = x_dist(rnd_engine);
  const Tile* t = this->resolve_coords(pt_x, pt_y, pt_k);
  assert(t != 0);
  return std::make_tuple( pt_x, pt_y, pt_k % this->layers.size() );
}

const Tile* Layer::get_rnd_tile(std::mt19937 rnd_engine) const {
std::uniform_int_distribution<std::size_t> t_dist(0, this->tiles.size() - 1);
  return this->tiles[t_dist(rnd_engine)];
}

std::size_t Layer::buffered_read(void* request_buf, std::size_t request_size, off_t request_offset) const {
  std::size_t pix_len = request_size / this->elem_size;
  std::size_t pix_off = request_offset / this->elem_size;
  std::size_t lxb = pix_off % this->dim.x;
  std::size_t lyb = pix_off / this->dim.x;
  std::size_t lxe = (pix_off + pix_len) % this->dim.x + 1;
  std::size_t lye = (pix_off + pix_len) / this->dim.x + 1;

  std::size_t buf_ptr = 0;                                  // number of bytes written to buffer
  std::size_t buf_pos = 0;                                  // current position in buffer
  char* rbuf = (char*) request_buf;

  /* layer read pattern:
   *    a   b       c   d
   *    0   3       12  16
   * 0 |    [***********|  lyb
   * 1 |****************|
   * 2 |***********]    |
   * 3 |                |  lye
   *   a = 0
   *   b = lxb
   *   c = lxe
   *   d = l->max.x
   */

//   std::cout << "READ: " << "x:" << lxb << "\ty:" << lyb << "\tx:" << lxe << "\ty:" << lye << "\tl:" << pix_len << std::endl;
//   std::cout << "READ: " << "lpx:" << this->pos.x << "\tlpy:" << this->pos.y << "\tlmx:" << this->max.x << "\tlmy:" << this->max.y << std::endl;
  bool multi_row = lyb != (lye - 1);
  std::size_t xmin, xmax, step = 16;
  for (std::size_t y = lyb; y < lye; y++) {
//     std::cout << "\trow: " << y << std::endl;
    xmin = 0;
    xmax = this->dim.x;
    if (y == lyb) {
      xmin = lxb;
    }
    if (y == (lye - 1)) {
      xmax = lxe;
    }
    
    
    /* tile read pattern:
     *    a   b       c   d
     *    0 2 3       12| 16
     * 0 |  | [*********|*|  lyb
     * 1 |**|***********|*|
     * 2 |**|********]  | |
     * 3 |  |           | |  lye
     *   a = 0
     *   b = lxb
     *   c = lxe
     *   d = l->max.x
     *    0    4    8    12   16   20   24
     *   [    |  *******|   *   * |   *     ]
     *   [    |    t1   |         |    t2   ]
     *                         posx   x
     */
//     std::cout << "\txmn: " << xmin << "\txmx: " << xmax << std::endl;
    for (std::size_t x = xmin; x < xmax; x += step) {
      const Tile* t = this->resolve_coords(x, y);
//       std::cout << "x: " << x << "\ty: " << y << "\tt: " << t << std::endl;
      if (t) {
//         std::cout << "\t\tx: " << x << std::endl;
        std::size_t txmin = 0,
                    txmax = xmax - t->pos.x,
                    ty = y - t->pos.y;
        if (x == xmin) {
          txmin = x - t->pos.x;
          txmax = xmax - t->pos.x;
        } else {
          ty = y - t->pos.y;
        }
        
        // set buffer position to the start of our offset into this tile
        buf_pos = ((y - lyb) * t->pos.x + txmin) * this->elem_size;
//         std::cout << buf_pos << std::endl;
        if (buf_pos >= pix_len) { break; }
        if (buf_ptr < buf_pos) {
//           std::cout << "nf:" << buf_pos - buf_ptr << std::endl;
          memset(&rbuf[buf_ptr], 0xff, buf_pos - buf_ptr);
          buf_ptr = buf_pos;
        }
        step = (t->max.x / 2) + 1;
        
        std::size_t read_len = txmax - txmin;
        if ((read_len + buf_pos/this->elem_size) > pix_len) {
          read_len = pix_len - (buf_pos / this->elem_size);
        }
        
//         std::cout << "bm:" << read_len + buf_pos/this->elem_size << " rl:" << read_len << std::endl;
        // read as much data from the tile as we can
        std::size_t np = t->read_row(&rbuf[buf_pos], read_len, xmin, y);
//         std::cout << np << std::endl;
//         assert (np <= read_len);
//         std::size_t np = 1024;
        buf_ptr += (np * this->elem_size);
      }
    }
  }
  
  buf_pos = pix_len * this->elem_size;
  if (buf_ptr < buf_pos) {
//     std::cout << "rp:" << buf_pos << " " << buf_ptr << " " <<  buf_pos - buf_ptr << std::endl;
    memset(&rbuf[buf_ptr], 0xff, buf_pos - buf_ptr);
    buf_ptr = buf_pos;
  }
  return buf_pos;
  /*
  std::size_t x = lxb;
  std::size_t step = 32;
  std::size_t xlim = 0;
  std::size_t xmin = lxb;
  // find the tiles (and indices) present in each row requested
  for (std::size_t y = lyb; y < lye && y < this->dim.y; y++) {
    if (y == lye) {
      xlim = lxe;
    } else {
      xlim = this->dim.x;
    }
     
     
    for (; x < xlim;) {
      const Tile* t = this->resolve_coords(x, y);
      if (t) {
        x += t->max.x;
        step = (t->max.x / 2) + 1;
        std::size_t ne = t->read((void*)(&((char*) request_buf)[buf_pos]), xlim - xmin, xmin, y);
        buf_pos += ne * this->elem_size;
        memset((void*)(&((char*) request_buf)[buf_pos]), 0xff, xlim - xmin - ne);
        buf_pos += (xlim - xmin - ne) * this->elem_size;
      } else {
        memset((void*)(&((char*) request_buf)[buf_pos]), 0xff, step);
        buf_pos += step * this->elem_size;
        x += step;
      }
    }
    x = 0; xmin = 0;
  }
  return buf_pos / this->elem_size;
  
  
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
  */
}

void swap_bytes(char* dbuf, char* sbuf, std::size_t buf_len, std::size_t type_len) {
  char* r_buf = sbuf + type_len, *l_buf = dbuf;
  while (l_buf < dbuf + buf_len) {
    for (int b = 0; b < type_len; b++) {
      *(l_buf++) = *(--r_buf);
    }
    r_buf += type_len * 2;
  }
}

std::size_t Tile::read_row(char* buf, std::size_t pix_len, std::size_t tx, std::size_t ty) const {
  std::size_t offset_beg = (ty * this->max.x + tx) * this->elem_size + map_start;
  std::size_t data_len = pix_len * this->elem_size;
  std::size_t offset_end = offset_beg + data_len;
  std::size_t row_end = (ty+1) * this->max.x * this->elem_size + map_start;

  if (offset_end > row_end) {
    data_len = row_end - offset_beg;
  }
  void* map_src = &((char*)this->map)[offset_beg];
  memcpy(buf, map_src, data_len);
  char tmp = 0;

  // swap endianness
  tmp = buf[0]; buf[0] = buf[3]; buf[3] = tmp;
  tmp = buf[1]; buf[1] = buf[2]; buf[2] = tmp;
  return data_len / this->elem_size;
}

// read len pixels into buf starting at (tx, ty)
// returns number of pixels read
std::size_t Tile::read_row_compressed(char* buf, std::size_t len, std::size_t tx, std::size_t ty) const {
  // Warning: compressed reads are not thread-safe
  // Note: the cfits API is indexed starting from 1
  int anynul = 0;
  int status = 0;
  long fpixel[2] = { (long)tx + 1, (long)ty + 1 };
  std::size_t req_len = len;
  if (fpixel[0] + req_len > this->dim.x) {
    req_len = this->max.x - fpixel[0];
  }
  fits_read_pix(this->fptr, TFLOAT, fpixel, req_len, 0, (float*)buf, &anynul, &status);
  return req_len;
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

} // namespace umap_fits_file_internal

namespace umap_fits_file {
// Public interface to umap_fits_file

void find_fits(std::string basedir, std::vector<std::string>& fitslist);

template <typename _pixel_type>
class umap_fits_cube {
public:

  using pixel_type = _pixel_type;
  using csv_map = std::unordered_map<std::string, std::string>;

  /// -------------------------------------------------------------------------------- ///
  /// Constructor
  /// -------------------------------------------------------------------------------- ///
  umap_fits_cube() = default;

  umap_fits_cube(const std::string& basepath, const char* datalist_file, const char* exposure_file) {
    size_t byte_per_element;
    m_cube.page_size = utility::umt_getpagesize();
    std::vector<std::string> fitslist;
    find_fits(basepath, fitslist);
    std::sort(fitslist.begin(), fitslist.end());

    csv_map csvdata;
    this->read_list_csv(csvdata, datalist_file);

    this->add_from_csv(csvdata, fitslist);
    this->size_x = m_cube.max.x + m_cube.pos.x;
    this->size_y = m_cube.max.y + m_cube.pos.y;
    this->size_k = m_cube.max.z + m_cube.pos.z;
    std::cout << m_timestamp_list.size() << " " << size_k << std::endl;
    assert(m_timestamp_list.size() == size_k);

    for (size_t i = 0; i < size_k; ++i) {
      m_exposuretime_list.push_back(40.0);
    }
//     this->cstore = new umap_fits_file_internal::CfitsStoreFile(&m_cube, m_cube.cube_size, m_cube.page_size);
//     const int prot = PROT_READ|PROT_WRITE;
//     int flags = UMAP_PRIVATE;
//     this->cstore->region = umap_ex(NULL, this->m_cube.cube_size, prot, flags, 0, 0, this->cstore);
//     if ( this->cstore->region == UMAP_FAILED ) {
//         ostringstream ss;
//         ss << "umap of " << this->m_cube.cube_size << " bytes failed for Cube";
//         perror(ss.str().c_str());
//         exit(1);
//     }
//     this->m_image_data = (pixel_type*)cstore->region;
  }

  ~umap_fits_cube() {
//     uunmap(this->cstore->region, this->m_cube.cube_size);
  }

  // prevent this object from being copied
  umap_fits_cube(const umap_fits_cube &) = delete;
  umap_fits_cube(umap_fits_cube &&) = delete;
  umap_fits_cube &operator=(const umap_fits_cube &) = delete;
  umap_fits_cube &operator=(umap_fits_cube &&) = delete;

  void add_from_csv(csv_map& csvdata, const std::vector<std::string>& fitslist) {
    double mjd_start;
    int count = 0;
    struct stat sbuf;
    using row_data_t = std::tuple<double, double, double, double, unsigned long, bool>;
    std::map<std::string, row_data_t> files_loaded;
    for (auto it = fitslist.begin(); it != fitslist.end(); it++) {
      std::string basename = it->substr(it->rfind("/") + 1);
      if ( stat(it->c_str(), &sbuf) == -1 ) {
        std::cerr << "Warning: file \"" << *it << "\" enumerated but does not exist." << std::endl;
        continue;
      }
      double psf = 2.0, mjd, ra = 0, dec = 0, noise = 100;
      unsigned long time = count * 4000;
      bool valid = false;
      row_data_t row_data = std::make_tuple(psf, ra, dec, noise, time, valid);
      if (!files_loaded.count(basename)) {
        auto _row = csvdata.find(basename);
        if (_row == csvdata.end()) {
          std::cerr << "Warning: file \"" << basename << "\" does not exist in data CSV." << std::endl;
          row_data = std::tie(psf, ra, dec, noise, time, valid);
        } else {
          std::string row(_row->second);
          if (sscanf(row.c_str(), "%lf,%lf,%lf,%lf,%lf", &psf, &mjd, &ra, &dec, &noise) == 5) {
            if (m_timestamp_list.size() == 0) {
              mjd_start = mjd;
            }
            time = std::round((mjd - mjd_start) * 24 * 60 * 60 * 100); // Hundreths of a second
            valid = true;
            row_data = std::tie(psf, ra, dec, noise, time, valid);
          } else {
            // NAN row
          }
        }
        files_loaded.emplace(basename, row_data);
        count += 1;
      }
      m_cube.add_tile( *it );
      files_loaded.emplace( basename, row_data );
    }
    for (auto kv : files_loaded) {
      double psf, ra, dec, noise; unsigned long time; bool valid = false;
//       std::tie(psf, ra, dec, noise, time, valid) = kv.second;
//       std::cout << kv.first << " " << psf << " " << ra << " " << dec << " " << noise << " " << time << " " << valid << std::endl;
      if (std::get<5>(kv.second)) {
        std::array<double, 2> ra_dec = {ra,  dec};
        m_timestamp_list.push_back(time);
        m_psf_list.push_back(psf);
        m_ra_dec_list.push_back(ra_dec);
        m_noise_list.push_back(noise);
      } else {
        // copy previous row data if NAN row
        m_timestamp_list.emplace_back(m_timestamp_list.back());
        m_psf_list.emplace_back(m_psf_list.back());
        m_ra_dec_list.emplace_back(m_ra_dec_list.back());
        m_noise_list.emplace_back(m_noise_list.back());
      }
    }
    m_cube.refit();
  }

  /// -------------------------------------------------------------------------------- ///
  /// Public methods
  /// -------------------------------------------------------------------------------- ///

  /// \brief Return frame size
  std::size_t frame_size() const {
    return m_cube.layer_size;
  }

  /// \brief Return cube size
  std::size_t cube_size() const {
    return m_cube.cube_size;
  }

  std::tuple<std::size_t, std::size_t, std::size_t> get_rnd_coord(std::mt19937 rnd_engine, double x_slope, double y_slope) const {
    auto intercept = this->m_cube.get_rnd_coord(rnd_engine);
    std::size_t k = std::get<2>(intercept);
    double time_offset = this->timestamp(k) - this->timestamp(0);
    std::size_t x = std::round(std::get<0>(intercept) - x_slope * time_offset);
    std::size_t y = std::round(std::get<1>(intercept) - y_slope * time_offset);
    return std::make_tuple(x, y, 0);
  }

  bool out_of_range(const std::size_t x, const std::size_t y, const std::size_t k) const {
    return (m_cube.resolve_coords(x, y, k) == 0);
  }

  /// \brief Returns the index of the given (x, y, k) in a 3D cube
  /// Returns -1 if the given xyz coordinate points at out side of the cube
  std::size_t index_in_cube(const ssize_t x, const ssize_t y, const ssize_t k) const {
    umap_fits_file_internal::Point real(x - m_cube.pos.x, y - m_cube.pos.y, k - m_cube.pos.z);
    if (real.x < 0 || size_x <= real.x || real.y < 0 || size_y <= real.y || real.z < 0 || size_k <= real.z) {
#if FITS_CUBE_VERBOSE_OUT_OF_RANGE
      std::cerr << "Cube index is out-of-range: "
              << "(" << real.x << ", " << real.y << ", " << real.z << ") is out of "
              << "(" << size_x << ", " << size_y << ", " << size_k << ")" << std::endl;
#endif
      return -1;
    }
    return m_cube.index(real);
  }

  /// \brief Returns a pixel value of the given x-y-k coordinate
  /// A returned value can be NaN value
  pixel_type get_pixel_value(const std::size_t x, const std::size_t y, const std::size_t k) const {
    const umap_fits_file_internal::Tile* t = m_cube.resolve_coords(x, y, k);
    if (!t) {
      return this->nan;
    }
    pixel_type rv = 0;
    t->read_row((char*)&rv, 1, x - t->pos.x, y - t->pos.y);
    return rv;
  }

  /// \brief Returns the size of cube (x, y, k) in tuple
  std::tuple<size_t, size_t, size_t> size() const {
    return std::make_tuple( size_x, size_y, size_k );
  }

  void import_hdu_data(const char* filename, std::vector<double>& dest, double* default_val) {
    if (filename != nullptr) {
      std::ifstream ifs(filename);
      if (!ifs.is_open()) {
        std::cerr << "Cannot open " << filename << std::endl;
        std::abort();
      }
      for (double timestamp; ifs >> timestamp;) {
        dest.emplace_back(timestamp);
      }
      if (dest.size() != this->size_k) {
        std::cerr << "#of lines in " << filename << " is not the same as #of fits files" << std::endl;
        std::abort();
      }
    } else {
      // If a list of timestamps is not given, assume that the difference between two frames is 1.0
      dest.resize(this->size_k);
      if (default_val != nullptr) {
        for (size_t i = 0; i < this->size_k; ++i) { dest[i] = *default_val; }
      } else {
        for (size_t i = 0; i < this->size_k; ++i) { dest[i] = i; }
      }
    }
  }

  // Function to read data info from a csv file
  // Reads timestamp, psf fwhm, (ra/dec), and background sky noise 
  void read_list_csv(csv_map& csvdata, const char* list_file_name) {
    if (list_file_name != nullptr) {
      std::ifstream ifs(list_file_name);
      if (!ifs.is_open()) {
        std::cerr << "Cannot open " << list_file_name << std::endl;
        std::abort();
      }
      std::string line;
      std::getline(ifs, line); //skip first row of header info
      char fn[32], fdata[256];
      while (std::getline(ifs, line)) {
        std::istringstream iss{ line };
        if (sscanf(iss.str().c_str(), "%[^,],%s", fn, fdata) == 2) {
          csvdata.emplace(std::string(fn), std::string(fdata));
        }
      }
    }
  }

  unsigned long timestamp(const size_t k) const {
    return m_timestamp_list[k];
  }

  double exposuretime(const size_t k) const {
    return m_exposuretime_list[k];
  }

  double psf(const size_t k) const {
    return m_psf_list[k];
  }

  std::array<double, 2> ra_dec(const size_t k) const {
    return m_ra_dec_list[k];
  }

  double noise(const size_t k) const {
    return m_noise_list[k];
  }

  /// -------------------------------------------------------------------------------- ///
  /// Public fields
  /// -------------------------------------------------------------------------------- ///

  std::size_t size_x;                                       // in pixels
  std::size_t size_y;                                       // in pixels
  std::size_t size_k;                                       // in pixels
  constexpr static const pixel_type nan = std::numeric_limits<double>::quiet_NaN();

private:
  /// -------------------------------------------------------------------------------- ///
  /// Private methods
  /// -------------------------------------------------------------------------------- ///

  /// -------------------------------------------------------------------------------- ///
  /// Private fields
  /// -------------------------------------------------------------------------------- ///

  pixel_type * m_image_data;
  std::vector<unsigned long> m_timestamp_list;    // an array of the timestamp of each frame.
  std::vector<double>        m_exposuretime_list; // an array of the exposure time of each image
  std::vector<double>        m_psf_list;          // an array of the psf fwhm of each image
  std::vector<std::array<double, 2>>     m_ra_dec_list;       // an array of ra/dec values for boresight of each image
  std::vector<double>        m_noise_list;        // an array of average background sky value (noise) for each image

  umap_fits_file_internal::Cube m_cube;
  umap_fits_file_internal::CfitsStoreFile* cstore;
};

/* Find all fits files in a base directory and all of its subdirectories */
void find_fits(std::string basedir, std::vector<std::string>& fitslist) {
  struct stat st;
  const char* basedir_c = basedir.c_str();
  if ((stat(basedir_c, &st) == -1) || !(st.st_mode & S_IFDIR)) {
    return;
  }

  DIR* dfd = opendir(basedir_c);
  struct dirent* entry = readdir(dfd);
  for (; entry; entry = readdir(dfd)) {
    std::string d_entry(entry->d_name);
    std::string d_path = basedir + "/" + d_entry;
    if (entry->d_name[0] == '.') { continue; }
    else if (d_entry.length() > 4 && d_entry.rfind(".fits") == (d_entry.size() - 5)) {
      fitslist.push_back(d_path);
    } else {
      find_fits(d_path, fitslist);
    }
  }
}

} // namespace umap_fits_file
} // namespace utility
#endif

