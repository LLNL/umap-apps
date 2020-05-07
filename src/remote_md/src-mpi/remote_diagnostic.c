/// An example client that monitors the atoms
/// and save analytics if required

#include "initAtoms.h"
#include "umap/umap.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
  if(argc!=3){
    printf("Usage: %s [nun_sim_procs] [max_time] \n\n", argv[0]);
    return 0;
  }

  int rank, num_procs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t umap_page_size = umapcfg_get_umap_page_size();
  int stride = umap_page_size/sizeof(real_t);
  int num_sim_procs = atoi(argv[1]);
  int num_time_steps = atoi(argv[2]);
  real_t* kinetic_energy = (real_t*) umap_network("kinetic_energy", NULL, num_sim_procs*umap_page_size);

  if ( kinetic_energy == UMAP_FAILED ) {
    printf("umap kinetic_energy from network-based datastore \n");
    return -1;
  }

  setbuf(stdout, NULL);
  int timestep = 0;
  printf("Monitoring Kinetic Energy ... \n");
  while (timestep<num_time_steps){

    printf("Time step = %d \n", timestep);
    int proc = 0;
    while( proc<num_sim_procs ){
      if(kinetic_energy[proc*stride+(timestep+1)]>0.0)
	{
	  printf("%f ", proc, kinetic_energy[proc*stride+(timestep+1)]);
	  proc ++;
	}
      else{
	umap_evict();
	sleep(5);
      }
    }
    printf("\n");
    timestep ++;
    
  }
  
  uunmap(kinetic_energy, 0);
  return 0;
}
