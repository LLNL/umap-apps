/// An example client that monitors the atoms
/// and save analytics if required

#include "initAtoms.h"
#include "umap/umap.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{

  int rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t umap_page_size = umapcfg_get_umap_page_size();
  real3* atoms_r = (real3*) umap_network("atoms_r", NULL, 332398592);

  if ( atoms_r == UMAP_FAILED ) {
    printf("umap atoms_r from network-based datastore \n");
    return -1;
  } 
  size_t numbytes = 0; //ds->get_size();
  int maxTotalAtoms = 332398592/sizeof(real3);


  srand(0);
  for (int i=0;i<10;i++){
    size_t pos = rand() % maxTotalAtoms;
    real3 tmp;
    tmp[0] = atoms_r[pos][0];
    tmp[1] = atoms_r[pos][1];
    tmp[2] = atoms_r[pos][2];
    printf("atoms_r[%zu]=%f %f %f \n", pos, tmp[0], tmp[1], tmp[2]);
    sleep(10);
  }
  
  uunmap(atoms_r, numbytes);
  return 0;
}
