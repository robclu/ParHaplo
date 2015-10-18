#include "tree_gpu.cpp"

TreeGpu::TreeGpu(binary_vector& data    , read_info_container& read_info, snp_info_container snp_info   ,   
                 const size_t   snps    , const size_t         reads    , const size_t       min_ubound ,
                 const size_t   device  )
: _data(data)            , _read_info(read_info) , _snp_info(snp_info)    , _haplotype(snps)  , 
  _alignments(reads)     , _tree(snps, reads)    , _snps(snps)            , _reads(reads)     ,
  _min_ubound(min_ubound), _device(device)
{
    // Move all the data to the device
    move_data();
}

~TreeGpu() 
{
    cudaFree(_tree.data);
    cudaFree(_tree.read_info);
    cudaFree(_tree.snp_info);
    cudaFree(_tree.nodes);
    cudaFree(_tree.haplotype);
    cudaFree(_tree.alignments);
    cudaFree(_tree.search_snps);
    cudaFree(_tree.aligned_reads);
    cudaFree(_selection_parameters);
}

void TreeGpu::move_data()
{
    cudaError_t error;
    
    // Copy data to device
    thrust::host_vector<small_type> data_gpu_format = _data.to_binary_vector();
    cudaMalloc((void**)&_tree.data, sizeof(small_type) * _data.size());
    cudaMemcpy(_tree.data, thrust::raw_pointer_cast(&_data_gpu_format[0]), 
                    sizeof(small_type) * data_gpu_format.size(), cudaMemcpyHostToDevice); 

    // Copy the read data to the device 
    cudaMalloc((void**)&_tree.read_info, sizeof(read_info_type) * _read_info.size());
    cudaMemcpy(_tree.read_info, thrust::raw_pointer_cast(&_read_info[0]),
                    sizeof(read_info_type) * _read_info.size(), cudaMemcpyHostToDevice); 
    
    // Copy SNP data to the device 
    cudaMalloc((void**)&_tree.snp_info, sizeof(snp_info_type) * _snp_info.size());
    cudaMemcpy(_tree.snp_info, thrust::raw_pointer_cast(&_snp_info[0]), 
                    sizeof(snp_info_type) * _snp_info.size(), cudaMemcpyHostToDevice); 
  
    // Allocate the haplotype and alignment data on the device 
    cudaMalloc((void**)&_tree.haplotype, sizeof(small_type) * _snps);
    cudaMalloc((void**)&_tree.alignments, sizeof(small_type) * _read_info.size());
    
    // Create a vector for the search snps
    thrust::host_vector<size_t> snp_indices;
    for (size_t i = 0; i < _snps; ++i) snp_indices.push_back(i);
    
    // Copy snp search info to the device
    cudaMalloc((void**)&_tree.search_snps, sizeof(size_t) * _snps);
    cudaMemcpy(_tree.search_snps, thrust::raw_pointer_cast(&snp_indices[0]), 
                sizeof(size_t) * snp_indices.size(), cudaMemcpyHostToDevice);
    
    // Create a vector for the aligned rows
    thrust::host_vector<size_t> aligned;
    for (size_t i = 0; i < _reads; ++i) aligned.push_back(i);
   
    // Copy aligned data to the device
    cudaMalloc((void**)&_tree.aligned_reads, sizeof(size_t) * _reads);
    cudaMemcpy(_tree.aligned_reads, thrust::raw_pointer_cast(&aligned[0]), 
                sizeof(size_t) * aligned.size(), cudaMemcpyHostToDevice);
    
    // Create data for bounds array
    cudaMalloc((void**)&_snp_bounds, sizeof(BoundsGpu) * _snps);

    // ------------------------------------ NODE ALLOCATION -------------------------------------------------
    
    // Allocate space on the device for the tree nodes
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
   
    // DEBUG 
    std::cout << "Free Memory Before : " << free_mem << "\n";
    
    // Determine the number of nodes to "safely allocate"
    size_t num_nodes = 0.7 * free_mem / sizeof(TreeNode);
    
    // Check that we can allocate the memory
    error = cudaMalloc((void**)&_tree.nodes, sizeof(TreeNode) * num_nodes);
    if (error != cudaSuccess) std::cout << "Error allocating memory\n";
    
    // DEBUG
    cudaMemGetInfo(&free_mem, &total_mem);
    std::cout << "Free Memory After: " << free_mem  << "\n";
    std::cout << "Num Nodes        : " << num_nodes << "\n";

    // ---------------------------------------- HEAP LIMIT ---------------------------------------------------
    
    const size_t heap_mb = 128;                         // Excessive, but just incase
    cudaThreadSetLimit(cudaLimitMallocHeapSize, heap_mb * 1024 * 1024);
    
    // --------------------------------------- TREE SEARCH ---------------------------------------------------

    
    // Invoke the kernel with just a single thread -- kernel spawns more threads
    search_tree<<<1,1>>>(_tree, _selection_parameters, _min_ubound, _device);
}

void TreeGpu::search_tree()
{
    // Number of snps that still need to be searched
    size_t unsearched_snps = snps - 1;
   
    // --------------------------------------- ROOT NODE ----------------------------------------------------
    
    // The first level of the tree (with the root node) is a little different
    // because it needs to do the setup, so invoke that kernel first
    map_root_node<<<1, 1>>>(_tree, _snp_bounds, _min_ubound, _device); 
    if (cudaSuccess != cudaGetLastError())  printf("Kernel Launch Error for Root Map\n"); 
    cudaDeviceSynchronize();
    
    // Now we can do the mapping of the unsearched nodes
    map_unsearched_snps<<<1, unsearched_snps>>>(_tree, _snp_bounds, 0);
    if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Unsearched SNP Map!\n");
    cudaDeviceSynchronize();
    __syncthreads();
    
    // And the reduction
   reduce_unsearched_snps<<<1, unsearched_snps>>>(_tree, unsearched_snps); 
    if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Unsearched SNP Reduce!\n");
    cudaDeviceSynchronize();
    __syncthreads();
    
    // ----------------------------------------- OTHER NODES ------------------------------------------------

    while () {
    }
    
    // Haplotype found
}


