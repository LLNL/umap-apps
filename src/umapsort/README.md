
# UMap Sort

## Map in a file of integers and then sort it

### Example

Sort an array of 96 GB stored in data_file using 4 threads.

```
drop_page_cache
free_mem 
env UMAP_PAGESIZE=$umap_psize ./umapsort -f ${SSD_MNT_PATH}/data_file -p $(((96*1024*1024*1024)/umap_psize)) -N 1 -t 4
```
