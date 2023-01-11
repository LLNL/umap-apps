// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#include "umap/umap.h"
#define MEM_MAP_WITH_UMAP 1
#define MEM_MAP_WITH_MMAP 2

_seq<char> memmapFromFileOffset(const char *filename, off_t offset, size_t length, void* addr, int mem_map_option)
{
  struct stat sb;
  int flags = O_RDWR; //O_RDONLY  O_RDWR
  if( mem_map_option == MEM_MAP_WITH_UMAP )
    flags |= O_DIRECT;

  int fd = open(filename, flags);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG (sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }

  if( addr && munmap(addr, length) == -1){
    printf("munmap to specified address at %p \n failed",addr);
    exit(-1);
  }

  char *p;
  if( mem_map_option == MEM_MAP_WITH_MMAP )
  {
    p = static_cast<char*>(mmap(addr, length, PROT_READ|PROT_WRITE, MAP_PRIVATE, fd, offset));
    if (p == MAP_FAILED) {
      perror("mmap failed");
      exit(-1);
    }
    cout << "mmapped at " << (void*)p << endl;
    if (close(fd) == -1) {
      perror("close failed");
      exit(-1);
    }
  }else{
    uint64_t umap_pagesize = umapcfg_get_umap_page_size();
    length = ((length-1)/umap_pagesize + 1 )*umap_pagesize;
    p = static_cast<char*>(umap(addr, length, PROT_READ, UMAP_PRIVATE, fd, offset)); //|PROT_WRITE
    if (p == UMAP_FAILED) {
      perror("umap failed");
      exit(-1);
    }
    //if(length<=4295229440UL) umap_fetch_and_pin(p, length);
    cout << "umapped at " << (void*)p << endl;
  }
  
  //profiler.mt_register_address((void*)p, length);

  //for Ligra, it needs to be mapped exactly to the requested virtual address
  if( addr && p!=addr ){
    printf("failed to be mapped to %p (%p)\n", addr, p);
    exit(-1);
  }
  
  return _seq<char>(p, length);
}

template <class vertex>
graph<vertex> readUncompressedGraph(char* iFile, bool symmetric, bool binary, bool mmap) {
  if(binary) 
    return readGraphFromBinary<vertex>(iFile,symmetric);
  else 
    return readGraphFromFile<vertex>(iFile,symmetric,mmap);
}

template <class vertex>
graph<vertex> generateCompressedSymmetricGraphStoreFromBinary(char* fname) 
{
  return NULL;
}

template <class vertex>
graph<vertex> generateCompressedSymmetricGraphStoreFromFile(char* fname) 
{
  return NULL;
}

template <class vertex>
graph<vertex> generateUncompressedSymmetricGraphStoreFromBinary(char* fname) 
{
  return NULL;
}

template <class symmetricVertex>
void generateUncompressedSymmetricGraphStoreFromFile(char* fname) 
{ 
  printf("generateUncompressedSymmetricGraphStoreFromFile:: start\n");
  string s(fname);
  std::string ss_path(s);

  words W;
  _seq<char> S = readStringFromFile(fname);
  W = stringToWords(S.A, S.n);
#ifndef WEIGHTED
  if (W.Strings[0] != (string) "AdjacencyGraph")
#else
  if (W.Strings[0] != (string) "WeightedAdjacencyGraph")
#endif
  {
    cout << "Bad input file" << endl;
    abort();
  }

  long len = W.m -1;
  long n = atol(W.Strings[1]);
  long m = atol(W.Strings[2]);
#ifndef WEIGHTED
  if (len != n + m + 2)
#else
  if (len != n + 2*m + 2)
#endif
  {
    cout << "Bad input file" << endl;
    abort();
  }else{
    printf("generateUncompressedSymmetricGraphStoreFromFile:: n = %ld m=%ld\n", n, m);
  }


  uintT* offsets = newA(uintT,n);
#ifndef WEIGHTED
  uintE* edges = newA(uintE,m);
#else
  intE* edges = newA(intE,2*m);
#endif


  {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
  {parallel_for(long i=0; i<m; i++) {
#ifndef WEIGHTED
      edges[i] = atol(W.Strings[i+n+3]);
#else
      edges[2*i] = atol(W.Strings[i+n+3]);
      edges[2*i+1] = atol(W.Strings[i+n+m+3]);
#endif
    }}
  //W.del(); // to deal with performance bug in malloc

  symmetricVertex* v = newA(symmetricVertex,n);

  {parallel_for (uintT i=0; i < n; i++) {
    uintT o = offsets[i];
    uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
    v[i].setOutDegree(l);
#ifndef WEIGHTED
    v[i].setOutNeighbors(edges+o);
#else
    v[i].setOutNeighbors(edges+2*o);
#endif
    }}

  free(offsets);

  /*
    Start CREATE_MMAP_FILE, the following saves graph edge structure into a file
  */
  std::stringstream ss_fname_edge;
  ss_fname_edge << ss_path << ".edge";
  const std::string str_edge = ss_fname_edge.str();
  const char*     fname_edge = str_edge.c_str();
  FILE * fp = fopen(fname_edge, "wb"); //write binary
  if (fp == NULL)  /* If an error occurs during the file creation */
  {
     cerr << "fopen() failed for " << fname_edge <<endl;
     exit(0);
  }

  fseek(fp, 0, SEEK_SET);
#ifndef WEIGHTED
  size_t elements_written = fwrite(edges, sizeof(uintE), m, fp); 
  if (elements_written != m )
    cerr << "fwrite(edges) failed" << endl;
  fseek(fp, sizeof(uintE)*m, SEEK_SET);
#else
  size_t elements_written = fwrite(edges, sizeof(uintE), m*2, fp); 
  if (elements_written != m*2 )
    cerr << "fwrite(weigthed edges) failed" << endl;
  fseek(fp, sizeof(uintE)*m*2, SEEK_SET);
#endif


  elements_written = fwrite(&m, sizeof(long), 1, fp); 
  if (elements_written != 1)
    cerr << "fwrite(m) failed "<<endl;

  //fseek(fp, sizeof(uintE)*m+sizeof(long), SEEK_SET);
  fseek(fp, 0L, SEEK_CUR);
  elements_written = fwrite(&n, sizeof(long), 1, fp); 
  if (elements_written != 1)
    cerr << "fwrite(n) failed "<<endl;

  //fseek(fp, sizeof(uintE)*m+sizeof(long)*2, SEEK_SET);
  fseek(fp, 0L, SEEK_CUR);
  uint64_t edges_ptr = (uint64_t) edges;
  elements_written = fwrite(&edges_ptr, sizeof(uint64_t), 1, fp); 
  if (elements_written != 1)
    cerr << "fwrite(edges_ptr) failed "<<endl;

  //truncate the file size
  //fseek(fp, sizeof(uintE)*m+sizeof(long)*2+sizeof(uint64_t), SEEK_SET);
  fseek(fp, 0L, SEEK_CUR);
  fputc('\0', fp);
  fclose(fp);

    
  /*
    Start saving graph vertex structure into a file
  */
  std::stringstream ss_fname_vertex;
  ss_fname_vertex << ss_path << ".vertex";
  const std::string str_vertex = ss_fname_vertex.str();
  const char*     fname_vertex = str_vertex.c_str();
  fp = fopen(fname_vertex, "wb"); //write binary
  if (fp == NULL)  /* If an error occurs during the file creation */
  {
     cerr << "fopen(vertex) failed for " << fname_vertex <<endl;
     exit(0);
  }
  fseek(fp, 0, SEEK_SET);
  elements_written = fwrite(v, sizeof(symmetricVertex), n, fp); 
  if (elements_written != n )
    cerr << "fwrite(vertex) failed" <<endl;
  
  //fseek(fp, sizeof(vertex)*n, SEEK_SET);
  fseek(fp, 0L, SEEK_CUR);
  uint64_t vertex_ptr = (uint64_t) v;
  elements_written = fwrite(&vertex_ptr, sizeof(uint64_t), 1, fp); 
  if (elements_written != 1)
    cerr << "fwrite(vertex_ptr) failed "<<endl;

  //truncate the file size
  //fseek(fp, sizeof(vertex)*n+sizeof(uint64_t), SEEK_SET);
  fseek(fp, 0L, SEEK_CUR);
  fputc('\0', fp);
  fclose(fp);

}

template <class vertex>
graph<vertex> memmapUncompressedSymmetricGraphStoreFromFile(char* fname, int mem_map_option) 
{

  string s(fname);
  std::string ss_path(s);

  // Graph Edge datastore
  std::stringstream ss_edge;
  ss_edge << ss_path << ".edge";
  const std::string str_edge = ss_edge.str();
  const char* fname_edge     = str_edge.c_str();

  /*
    Start of reading the graph metadata from edge file/datastore : num_edges, num_vertices, edge_ptr
  */
  FILE *fp = fopen(fname_edge, "rb");
  int fd = open(fname_edge, O_RDONLY);
  struct stat sb;
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }

  off_t offset = sb.st_size - 1 - (sizeof(long)*2+sizeof(uint64_t));
  //cout << "edge file_size " << sb.st_size << " offset " << offset << endl;

  long buffer[2];
  fseek(fp, offset, SEEK_SET);
  fread(buffer, sizeof(long), 2, fp);
  size_t num_edges    = buffer[0];
  size_t num_vertices = buffer[1];

  // check the file size and the number of edges are consistent
  if (offset != num_edges*sizeof(uintE)) {
    cerr << "edge file validation failed, compiled with wrong options?" << endl
         << "file offset " << offset << ", expected " << num_edges*sizeof(uintE) << endl;
    exit(1);
  }

  // read in the saved edge pointer
  fseek(fp, offset+sizeof(long)*2, SEEK_SET);
  uint64_t edge_ptr;
  fread(&edge_ptr, sizeof(uint64_t), 1, fp);
  void* p = (void*) edge_ptr;
  fclose(fp);

  /*
    Start of reading the graph metadata from vertex file/datastore : vertex_ptr
  */
  std::stringstream ss_vertex;
  ss_vertex << ss_path << ".vertex";
  const std::string str_vertex = ss_vertex.str();
  const char*  fname_vertex    = str_vertex.c_str();

  fp = fopen(fname_vertex, "rb");
  fd = open(fname_vertex, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
    
  // check the file size and the number of vertices are consistent
  size_t total_vertex_bytes = sb.st_size - 1 - sizeof(uint64_t);
  size_t element_bytes = sizeof(symmetricVertex);
  if ( total_vertex_bytes % element_bytes || total_vertex_bytes/element_bytes != num_vertices ) {
    cerr << "vertex file validation failed" << endl
         << "total_vertex_bytes " << total_vertex_bytes << ", expected " << num_vertices*element_bytes << endl;
    exit(1);
  }

  // read in the saved vertex pointer
  fseek(fp, sb.st_size -1 -sizeof(uint64_t), SEEK_SET);
  uint64_t vertex_ptr;
  fread(&vertex_ptr, sizeof(uint64_t), 1, fp);
  void* v_p = (void*) vertex_ptr;

  fclose(fp);

  /*
    Start of memory-mapping vertex and edge datastores
  */
  size_t edge_region_size = sizeof(uintE) * num_edges;
  off_t  edge_file_offset = 0;

  _seq<char> seq_edge = memmapFromFileOffset(fname_edge, edge_file_offset, edge_region_size, p, mem_map_option);
  uintE *edge_region  = (uintE*) seq_edge.A;
  if (edge_region == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  printf("Mem-mapped edges to %p of %lu bytes \n", edge_region, edge_region_size);


  size_t vertice_region_size = sizeof(symmetricVertex) * num_vertices;
  size_t vertice_file_offset = 0;
  _seq<char> seq_vertex = memmapFromFileOffset(fname_vertex, vertice_file_offset, vertice_region_size, NULL, mem_map_option);//vertex does not need to be mapped to the same address
  symmetricVertex* vertex_region = (symmetricVertex*) seq_vertex.A;
  if (vertex_region == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  printf("Mem-mapped vertex to %p of %lu bytes \n", vertex_region, vertice_region_size);

  long n = num_vertices;
  long m = num_edges;
  uintE* edges = (uintE*) edge_region;
  symmetricVertex* v = (symmetricVertex*) vertex_region;

  Uncompressed_Mem<vertex>* mem = new Uncompressed_Mem<vertex>(v,n,m,edges);
  cout << "Mem-map Uncompressed Symmetric UNWEIGHTED AdjacencyGraph, edges(m)=" << m <<", vertices(n)=" << n << endl;

  return graph<vertex>(v,n,m,mem);

}
