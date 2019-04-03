/*
 * This file is part of UMAP.  For copyright information see the COPYRIGHT file
 * in the top level directory, or at
 * https://github.com/LLNL/umap/blob/master/COPYRIGHT. This program is free
 * software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.  This program is distributed in
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
 * the terms and conditions of the GNU Lesser General Public License for more
 * details.  You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef _COMMANDLINE_HPP
#define _COMMANDLINE_HPP

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif // _GNU_SOURCE

#include <stdint.h>
#include <iostream>     // cout/cerr
#include <list>
#include <unistd.h>     // getopt()
#include <getopt.h>     // duh...
#include "umap/umap.h"
#include <omp.h>

namespace utility {
typedef struct {
  int initonly;       // Just perform initialization, then quit
  int noinit;         // Init already done, so skip it
  int usemmap;
  int shuffle;

  long pagesize;
  uint64_t numpages;
  uint64_t numthreads;
  uint64_t bufsize;
  uint64_t uffdthreads;
  uint64_t pages_to_access;  // 0 (default) - access all pages
  char const* filename; // file name or basename
  char const* dirname; // dir name or basename
} umt_optstruct_t;

static char const* DIRNAME = "/mnt/intel/";
static char const* FILENAME = "abc";

static const uint64_t NUMPAGES = 10000000;
static const uint64_t NUMTHREADS = 2;
static const uint64_t BUFFERSIZE = 16;

static std::list< void(*)(char*) > umt_usageprinters;

/* If this function returns a non-zero condition, then print the usage function and quit */
static void default_usage(char* pname)
{
#ifdef _OPENMP
  uint64_t num_threads = omp_get_num_threads();
#else
  uint64_t num_threads = NUMTHREADS;
#endif
  std::cerr
  << "Usage: " << pname << " [--initonly] [--noinit] [--directio]"
  <<                       " [--usemmap] [-p #] [-t #] [-b #] [-f name]\n\n"
  << " --help                 - This message\n"
  << " --initonly             - Initialize file, then stop\n"
  << " --noinit               - Use previously initialized file\n"
  << " --usemmap              - Use mmap instead of umap\n"
  << " --shuffle              - Shuffle memory accesses (instead of sequential access)\n"
  << " -p # of pages          - default: " << NUMPAGES << std::endl
  << " -t # of threads        - default: " << num_threads << std::endl
  << " -u # of uffd threads   - default: " << umap_cfg_get_uffdthreads() << " worker threads\n"
  << " -b # page buffer size  - default: " << umap_cfg_get_bufsize() << " Pages\n"
  << " -a # pages to access   - default: 0 - access all pages\n"
  << " -f [file name]         - backing file name.  Or file basename if multiple files\n"
  << " -d [directory name]    - backing directory name.  Or dir basename if multiple dirs\n"
  << " -P # page size         - default: " << umap_cfg_get_pagesize() << std::endl;
}

static int default_getoptions(void* optstruct, int argc, char *argv[])
{
  utility::umt_optstruct_t* testops = (utility::umt_optstruct_t*) optstruct;
  int c;
  char *pname = argv[0];

  testops->initonly = 0;
  testops->noinit = 0;
  testops->usemmap = 0;
  testops->shuffle = 0;
  testops->pages_to_access = 0;
  testops->numpages = NUMPAGES;
  testops->numthreads = NUMTHREADS;
  testops->bufsize = umap_cfg_get_bufsize();
  testops->uffdthreads = umap_cfg_get_uffdthreads();
  testops->filename = FILENAME;
  testops->dirname = DIRNAME;
  testops->pagesize = umap_cfg_get_pagesize();

#ifdef _OPENMP
  testops->numthreads = omp_get_num_threads();
#else
  testops->numthreads = NUMTHREADS;
#endif

  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"initonly",  no_argument,  &testops->initonly, 1 },
      {"noinit",    no_argument,  &testops->noinit,   1 },
      {"usemmap",   no_argument,  &testops->usemmap,  1 },
      {"shuffle",   no_argument,  &testops->shuffle,  1 },
      {"help",      no_argument,  NULL,  0 },
      {0,           0,            0,     0 }
    };

    c = getopt_long(argc, argv, "p:t:f:b:d:u:a:P:", long_options, &option_index);
    if (c == -1)
      break;

    switch(c) {
      case 0:
        if (long_options[option_index].flag != 0)
          break;
        return -1;

      case 'P':
        if ((testops->pagesize = strtol(optarg, nullptr, 0)) > 0) {
          if (umap_cfg_set_pagesize(testops->pagesize) < 0) {
            return -1;
          }
          break;
        }
        return -1;
      case 'p':
        if ((testops->numpages = strtoull(optarg, nullptr, 0)) > 0)
          break;
        return -1;
      case 't':
        if ((testops->numthreads = strtoull(optarg, nullptr, 0)) > 0)
          break;
        else return -1;
      case 'b':
        if ((testops->bufsize = strtoull(optarg, nullptr, 0)) > 0)
          break;
        else return -1;
      case 'u':
        if ((testops->uffdthreads = strtoull(optarg, nullptr, 0)) > 0)
          break;
        else return -1;
      case 'a':
        testops->pages_to_access = strtoull(optarg, nullptr, 0);
        break;
      case 'd':
        testops->dirname = optarg;
        break;
      case 'f':
        testops->filename = optarg;
        break;
    }
  }

  if (testops->numpages < testops->pages_to_access) {
    std::cerr << "Invalid -a argument " << testops->pages_to_access << "\n";
    return -1;
  }

  /*
   * Note: Care must be taken when configuring the number of threads
   * and the buffer size of umap.  When the buffer size is set, it
   * apportions the buffer evenly to the umap threads.  So setting the
   * buffer size requires that the number of threads be set properly
   * first.
   */
  if (testops->uffdthreads != umap_cfg_get_uffdthreads())
    umap_cfg_set_uffdthreads(testops->uffdthreads);

  #ifdef _OPENMP
    omp_set_num_threads(testops->numthreads);
  #endif
  umap_cfg_set_bufsize(testops->bufsize);

  return 0;
}

static int umt_opts_argc;
static char** umt_opts_argv;

void umt_handle_options(void* optstruct, int(*handler)(void*, int, char*[]), void(*usage)(char*) );

void umt_getoptions(utility::umt_optstruct_t* testopts, int argc, char *argv[])
{
  umt_opts_argc = argc;
  umt_opts_argv = argv;
  umt_usageprinters.push_front(&default_usage);
  umt_handle_options((void*) testopts, &default_getoptions, (void(*)(char*)) 0);
}

void umt_handle_options(void* optstruct, int(*handler)(void*, int, char*[]), void(*usage)(char*) )
{
  optind = 1;
  opterr = 0;
  if (usage) {
    umt_usageprinters.push_back(usage);
  }
  if (handler(optstruct, umt_opts_argc, umt_opts_argv)) {
    for (auto it = umt_usageprinters.begin(); it != umt_usageprinters.end(); it++) {
      (*it)(umt_opts_argv[0]);
      std::cout << std::endl;
    }
  }
}

long umt_getpagesize(void)
{
  return umap_cfg_get_pagesize();
}
}
#endif // _COMMANDLINE_HPP
