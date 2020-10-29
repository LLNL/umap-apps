#!/bin/bash

# -------------------------------------------------------- #
# Usage
# -------------------------------------------------------- #
# cd path/to/build dir/src/bfs/
# sh path/to/umap-apps/src/bfs/run_bfs_bench.sh

# -------------------------------------------------------- #
# Configuration
# -------------------------------------------------------- #
scale=18
out_file_prefix="bfs_s${scale}"
device_name="pmem0"
graph_file_path="/mnt/pmem/pm0/graph_s${scale}"
# -------------------------------------------------------- #

# -------------------------------------------------------- #
# Functions
# -------------------------------------------------------- #
set_up_variables() {
    num_vertices=$((2**${scale}))
    num_edges=$((${num_vertices}*32))

    base_bench_options="-n${num_vertices} -m${num_edges} -g${graph_file_path}"

    default_num_ra_blocks=$(($(cat "/sys/block/${device_name}/queue/read_ahead_kb")*2))
    echo "Default num readahead blocks = ${default_num_ra_blocks}"
}

execute_command() {
    echo "$@" |& tee -a ${out_file}
    time "$@" |& tee -a ${out_file}
}

used_gcc_version() {
    ret=$(strings $1 | grep "GCC")
    echo ${ret}
}

run() {
  # ---- Generate output file name ---- #
  if [ $usemmap -eq 1 ]; then
    ADDITIONAL_OPTION=" -s " # use system mmap
    out_file="${out_file_prefix}_m${usemmap}_t${num_app_threads}_ra${read_ahead_size}.log"
  else
    ADDITIONAL_OPTION=""
    out_file="${out_file_prefix}_m${usemmap}_t${num_app_threads}_f${umap_page_fillers}_e${umap_page_evictors}_h${umap_high_evict}_l${umap_low_evict}_p${umap_page_size}_r${umap_read_ahead}.log"
  fi
  date | tee ${out_file}

  # ---- Print some system information ---- #
  local source_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
  ${source_dir}/../tools/pretest_config >> ${out_file}
  df -lh >> ${out_file}
  local gcc_version=$(used_gcc_version ./run_bfs)
  echo "USED GCC is " ${gcc_version} >> ${out_file}
  echo "" | tee -a ${out_file}

  # ---- Set some environmental variables ---- #
  if [ $usemmap -eq 0 ]; then
      export UMAP_PAGE_FILLERS=$umap_page_fillers
      export UMAP_PAGE_EVICTORS=$umap_page_evictors
      export UMAP_EVICT_HIGH_WATER_THRESHOLD=$umap_high_evict
      export UMAP_EVICT_LOW_WATER_THRESHOLD=$umap_low_evict
      export UMAP_PAGESIZE=$umap_page_size
      export UMAP_READ_AHEAD=$umap_read_ahead
      env | grep "UMAP" |& tee -a ${out_file}
  fi
  export OMP_NUM_THREADS=$num_app_threads
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}" |& tee -a ${out_file}

  # ---- Run the benchmark ---- #
  execute_command ./run_bfs ${base_bench_options} ${ADDITIONAL_OPTION}
  echo "" |& tee -a ${out_file}

  date | tee -a ${out_file}
}


# -------------------------------------------------------- #
# Run benchmark varying configuration
# -------------------------------------------------------- #
main() {
    set_up_variables

    # ---- Run benchmark with mmap ---- #
    usemmap=1
    for num_app_threads in 48 96; do
        for read_ahead_size in 0 256; do
            /home/perma/change_readahead ${read_ahead_size} /dev/${device_name}
            run
            /home/perma/change_readahead ${default_num_ra_blocks} /dev/${device_name}
        done
    done

    # ---- Run benchmark with umap ---- #
    K=1024
    M=$((K*K))
    G=$((K*K*K))

    usemmap=0

    for num_app_threads in 48 96 #num_app_threads
    do
      for umap_page_fillers in 32 #umap_page_fillers
      do
        for umap_page_evictors in 16 #umap_page_evictors
        do
          for umap_high_evict in 90 #umap_high_evict
          do
            for umap_low_evict in 70 #umap_low_evict
            do
              for umap_page_size in $((4*K)) $((16*K)) $((64*K)) $((256*K)) $((1*M)) $((4*M)) #umap_page_size
              do
                for umap_read_ahead in 0 #umap_read_ahead
                do
                  run
                done
              done
            done
          done
        done
      done
    done
}

main "$@"