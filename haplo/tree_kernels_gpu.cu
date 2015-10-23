#include "cuda.h"
#include "debug.h"
#include "math.h"
#include "tree_internal.h"

namespace haplo {

struct BoundsGpu {

    size_t lower;               // The lower value for the bound    
    size_t upper;               // The upper value for the bound
    size_t diff;                // The difference between the upper and lower bound
    size_t index;               // The snp index the bound represents
    size_t offset;              // The offset in the array of the bound
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Overload for the assingment operator 
    // ------------------------------------------------------------------------------------------------------
    __device__
    void operator=(BoundsGpu& other) 
    {
        lower           = other.lower;
        upper           = other.upper;
        diff            = other.diff;
        index           = other.index;
        offset          = other.offset;
    }
};


// Maps all the unsearched snps to an array of BoundsGpu structs which can then be reduces
__global__ 
void map_unsearched_snps(internal::Tree tree, BoundsGpu* snp_bounds, const size_t last_searched_snp,
                         const size_t start_node_idx)
{
    const size_t    thread_id     = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t    ref_haplo_idx = tree.search_snps[thread_id + last_searched_snp + 1];
    const TreeNode* comp_node     = tree.node_ptr(start_node_idx);
    size_t same = 0, opp = 0;
    
    // Index of the haplotype for the comparison
    snp_bounds[thread_id].index  = ref_haplo_idx;
    snp_bounds[thread_id].offset = thread_id + last_searched_snp + 1;

    // Go through all the set nodes    
    for (size_t i = last_searched_snp + 1; i >= 1; --i) {
        const size_t alignments       = tree.alignment_offsets[comp_node->haplo_idx + tree.reads];
        const size_t alignment_offset = tree.alignment_offsets[comp_node->haplo_idx] - 1;        
 
        for (size_t j = alignment_offset; j < alignment_offset + alignments; ++j) {
            size_t row_offset = tree.read_info[tree.aligned_reads[j]].offset();
            size_t read_start = tree.read_info[tree.aligned_reads[j]].start_index();
            size_t read_end   = tree.read_info[tree.aligned_reads[j]].end_index();
#ifdef DEBUG 
            if (thread_id == 0) {
    printf("A : %i\n", tree.aligned_reads[j]);
            }
#endif            
            // If the reads cross the snp sites
            if (read_start <= ref_haplo_idx && read_start <= comp_node->haplo_idx &&
                read_end   >= ref_haplo_idx && read_end   >= comp_node->haplo_idx ) {
                
                // Values of the two elements
                uint8_t ref_value  = tree.data[row_offset + (ref_haplo_idx - read_start)];
                uint8_t comp_value = tree.data[row_offset + (comp_node->haplo_idx - read_start)]; 
             
                if (comp_value == ref_value && ref_value <= 1) ++same;
                else if (comp_value != ref_value && ref_value <= 1 && comp_value <= 1) ++opp;
            }
        }
        if (i > 1) {
            // Go back up the tree
            comp_node = tree.node_ptr(comp_node->root_idx);
        }
    }
    
    // Format the bounds 
    snp_bounds[thread_id].lower = min((unsigned int)same, (unsigned int)opp);
    snp_bounds[thread_id].upper = max((unsigned int)same, (unsigned int)opp);
    snp_bounds[thread_id].diff  = snp_bounds[thread_id].upper - snp_bounds[thread_id].lower;
    __syncthreads();

#ifdef DEBUG   
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        printf("MUS : UNSEARCHED SNPS : %i\n", last_searched_snp);
        for (size_t i = 0; i < tree.snps; ++i) 
            printf("%i ", snp_bounds[i].diff);
        printf("\n");
        for (size_t i = 0; i < tree.snps; ++i) 
            printf("%i ", snp_bounds[i].index);
        printf("\n");
    }
#endif
}

// Checks which of the two snps is more "vaulable", first by the bounds diff
// parameters, and then by the snp index
__device__
bool more_valuable(BoundsGpu* snp_one, BoundsGpu* snp_two)
{
    // First check which has the greater differenece between upper and lower
    return snp_one->diff > snp_two->diff 
                         ? true 
                         : snp_two->diff > snp_one->diff 
                            ? false 
                            : snp_one->index < snp_two->index 
                                ? true : false;
}

__device__ 
void swap_search_snp_indices(internal::Tree& tree, size_t last_searched_snp, const BoundsGpu* const bound)
{
#ifdef DEBUG
    printf("SSSI : NEXT UNSEARCHED : %i\n", last_searched_snp + 1);
    printf("SSSI : SWAP INDEX      : %i\n", bound->index);
#endif 
    ++last_searched_snp;
    // If the value to swap is not already the value
    if (bound->index != tree.search_snps[last_searched_snp]) {
        const size_t temp                   = tree.search_snps[last_searched_snp];
        tree.search_snps[last_searched_snp] = bound->index;
        tree.search_snps[bound->offset]     = temp;
    }
}

// "Reduce" function for the list of selection params to determine the best param 
//  The reduction is done such that the resulting bound has the highest diff and 
//  lowest index (this is to mean that the snp is most correlated to the snps already
//  searched and is also closest to them (hence the lowest index))
__global__ 
void reduce_unsearched_snps(internal::Tree tree, BoundsGpu* snp_bounds, const size_t last_searched_snp,
                            const size_t elements)
{
    const size_t thread_id      = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t reductions     = static_cast<size_t>(ceil(log2(static_cast<double>(elements))));
    size_t reduction_threads    = elements - 1, snp_idx_other = 0;
    
    while (reduction_threads > 1) {
        // Only the first half of the threads do the reduction 
        snp_idx_other = thread_id + (reduction_threads / 2);
        
        if (thread_id < (reduction_threads / 2)) {
            // If the more rightward bound is more valuable 
            if (!more_valuable(&snp_bounds[thread_id], &snp_bounds[snp_idx_other])) {
                // We need ro replace the left value withe the right one
                BoundsGpu temp            = snp_bounds[thread_id];
                snp_bounds[thread_id]     = snp_bounds[snp_idx_other];
                snp_bounds[snp_idx_other] = temp;
            }
        }
        // If there were an odd number of elements, the last one just fetches and moves a value 
        if (reduction_threads % 2 == 1) {
            if (thread_id == (reduction_threads / 2)) {
                // There is an odd number of elements in the array,
                // The last thread just needs to move a value 
                BoundsGpu temp            = snp_bounds[thread_id];
                snp_bounds[thread_id]     = snp_bounds[snp_idx_other];
                snp_bounds[snp_idx_other] = temp;
            }
            reduction_threads /= 2; reduction_threads += 1;
        } else reduction_threads /= 2;
        __syncthreads();
    }
    // The first thread needs to swap the search result
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        swap_search_snp_indices(tree, last_searched_snp, snp_bounds); 
    }
    __syncthreads();
        
#ifdef DEBUG
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        printf("RUS : SIZEOF BOUNF    : %i\n", sizeof(snp_bounds[0]));
        printf("RUS : UNSEARCHED SNPS : %i\n", elements);
        for (size_t i = 0; i < tree.snps; ++i) 
            printf("%i ", snp_bounds[i].diff);
        printf("\n");
        for (size_t i = 0; i < tree.snps; ++i) 
            printf("%i ", snp_bounds[i].index);
        printf("\n");
        for (size_t i = 0; i < tree.snps; ++i) 
            printf("%i ", tree.search_snps[i]);
        printf("\n");
    }
#endif
}

// Swaps two nodes in the node array
__device__ 
void swap(TreeNode* left, TreeNode* right) 
{
    TreeNode temp = *left;
    *left          = *right;
    *right         = temp;
}

// Updates two nodes
__device__
void update(TreeNode* left, TreeNode* right) 
{
    // Check if the node with the higher upper bound needs to be pruned
    if (left->min_ubound <= right->lbound) {
        const size_t right_pruned = right->pruned;
        left->pruned             += right_pruned;
        right->pruned            -= right_pruned;
        
        // If the right node has not been pruned yet 
        if (right->prune == 0) {
            ++left->pruned;
            right->prune = 1;
        }
        
#ifdef DEBUG 
        printf("U : LPRUNED : %i\n", left->pruned);
        printf("U : RPRUNED : %i\n", right->pruned);
#endif
    }
    // Update the minimun upper bounds
    left->min_ubound  = static_cast<size_t>(min(static_cast<double>(left->min_ubound) ,  
                                                static_cast<double>(right->min_ubound)));
    right->min_ubound = static_cast<size_t>(min(static_cast<double>(left->min_ubound) , 
                                                static_cast<double>(right->min_ubound)));
}

// Compares two tree nodes 

// "Reduces" the nodes, removing all the nodes with a lower bounds > than the highest upper bound
// by movind them to the end of the array and then returning the index of the end of the array
// also sorts the nodes by lower bound -- requires 2*log_2(N) operations

// Sorts a sub array bitonically 
__device__ 
void bitonic_out_in_sort(internal::Tree& tree       , const size_t start_offset    ,  
                         const size_t    block_index, const size_t num_nodes       ,
                         const size_t    total_nodes                               )
{
    const size_t start_id     = start_offset + block_index * num_nodes;
    const size_t thread_id    = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t node_id      = thread_id % num_nodes;
    const size_t comp_node_id = num_nodes - node_id - 1;
    
    if ((node_id < num_nodes / 2) && (block_index * num_nodes + comp_node_id < total_nodes)) {
        // First node compares with last node, second with second last ...
        if (tree.node_ptr(start_id  + comp_node_id)->lbound  < 
            tree.node_ptr(start_id + node_id)->lbound        ) {
                // Then the nodes need to be swapped
                swap(tree.node_ptr(start_id + comp_node_id), 
                     tree.node_ptr(start_id + node_id)    );
        }
    }
}


__device__
void bitonic_out_out_sort(internal::Tree& tree       , const size_t start_offset  ,
                          const size_t    block_index, const size_t num_nodes     ,
                          const size_t    total_nodes                             )
{
    const size_t start_id     = start_offset + block_index * num_nodes;
    const size_t thread_id    = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t node_id      = thread_id % num_nodes;
    const size_t comp_node_id = node_id + (num_nodes / 2);
    
    if ((node_id < num_nodes / 2) && (block_index * num_nodes + comp_node_id < total_nodes)) {
        // First node compares with the node __num_nodes__ ahead of it
        if (tree.node_ptr(start_id + comp_node_id)->lbound < 
            tree.node_ptr(start_id + node_id)->lbound      ) {
                // Then the nodes need to be swapped
                swap(tree.node_ptr(start_id + comp_node_id), 
                     tree.node_ptr(start_id + node_id)     );
        }
    }    
}

// Reduces a level using a parallel bitonic sort
__global__
void reduce_level(internal::Tree tree, const size_t start_node_idx, const size_t num_nodes)
{
    const size_t passes         = static_cast<size_t>(ceil(log2(static_cast<double>(num_nodes))));
    const size_t thread_id      = blockIdx.x * blockDim.x + threadIdx.x;
    size_t       block_size     = 2;

    for (size_t pass = 0; pass < passes; ++pass) {
        // We need a logarithmically decreasing number of out-in passes 
        if (thread_id % block_size < block_size / 2) 
            bitonic_out_in_sort(tree, start_node_idx, thread_id / block_size, block_size, num_nodes);
        __syncthreads();
        
        // Now we need pass number of out-out bitonic sorts
        size_t out_out_block_size = block_size / 2; 
        for (size_t i = 0; i < pass; ++i) {
            if (thread_id % out_out_block_size < out_out_block_size / 2) {
                bitonic_out_out_sort(tree, start_node_idx, thread_id / out_out_block_size, 
                                     out_out_block_size , num_nodes );
            }
            out_out_block_size /= 2;
            __syncthreads();
        }
        block_size *= 2;
    }
#ifdef DEBUG
            if (thread_id == 0) {
                printf("A : \n");
                for (size_t j = start_node_idx; j < start_node_idx + num_nodes; ++j) {
                    printf("%i ", tree.node_ptr(j)->lbound);
                } printf("\n"); 
            }
#endif   
}

// Updates the values of the alignments for a specific node
__device__
void add_node_alignments(internal::Tree& tree       , const size_t* const last_unaligned_idx, 
                         const size_t    node_idx   , const size_t* const alignment_offset  )
{
    TreeNode* const node            = tree.node_ptr(node_idx);
    const size_t thread_id          = blockIdx.x * blockDim.x + threadIdx.x;
    size_t alignments               = tree.alignment_offsets[node->haplo_idx + tree.reads];
    size_t alignment_start          = tree.alignment_offsets[node->haplo_idx];
    const size_t offset             = alignments * thread_id + *alignment_offset;
    size_t read_offset, elem_pos; uint8_t element;
    
    // If there is a valid index for the alignment start
    if (alignment_start != 0) {
        --alignment_start;
    
        // Set the offset in the read_values array of the start of the reads for this node
        node->align_idx = offset;
        
        // For all the alignments for the snp
        for (size_t i = alignment_start; i < alignment_start + alignments; ++i) {
            // Find the element in the data array
            read_offset = tree.read_info[tree.aligned_reads[i]].offset();
            elem_pos    = node->haplo_idx - tree.read_info[tree.aligned_reads[i]].start_index();
            element     = tree.data[read_offset + elem_pos];
            
            // Determine the alignment values of the node 
            if ((element == 0 && node->value == 0) || (element == 1 && node->value == 1)) {
                tree.read_values[offset + i - alignment_start] = 1;
            } else if ((element == 0 && node->value == 1) || (element == 1 && node->value == 0)) {
                tree.read_values[offset + i - alignment_start] = 0;
            }
        }
    }
}

// Updates the indices which are have been aligned so far
__device__ 
void update_global_alignments(internal::Tree& tree, size_t* last_unaligned_idx, const size_t haplo_idx)   
{
    __shared__ size_t  found   ; __shared__ size_t read_offset; __shared__ size_t  elem_pos     ;
    __shared__ size_t  elements; __shared__ size_t last_index ; __shared__ uint8_t element_value;
    
    elements   = tree.snp_info[haplo_idx].elements();
    last_index = *last_unaligned_idx;
    found      = 0; 
    
    // Go through all the unaligned reads
    for (size_t i = last_index; i < tree.reads && found < elements; ++i) {
        if (tree.aligned_reads[i] >= tree.snp_info[haplo_idx].start_index() &&
            tree.aligned_reads[i] <= tree.snp_info[haplo_idx].end_index()   ) {
            // If the read crosses the snp site
            if (tree.read_info[tree.aligned_reads[i]].start_index() <= haplo_idx &&
                tree.read_info[tree.aligned_reads[i]].end_index()   >= haplo_idx ) { 
                
                // Get the offset in memory of the start of the read
                read_offset     = tree.read_info[tree.aligned_reads[i]].offset();
                elem_pos        = haplo_idx - tree.read_info[tree.aligned_reads[i]].start_index();
                element_value   = tree.data[read_offset + elem_pos];
    
                // If the element is valid then add the index of the alignment by swapping 
                if (element_value <= 1) {
                    size_t temp_value               = tree.aligned_reads[i];
                    tree.aligned_reads[i]           = tree.aligned_reads[last_index];
                    tree.aligned_reads[last_index]  = temp_value;
                    ++last_index; ++found;
                }
            }
        }
    }
    // Update the last index 
    *last_unaligned_idx = last_index;
}

// Does the mapping for the leaf node -- calculates the alignments
__device__ 
void map_leaf_bounds(internal::Tree& tree, const size_t last_searched_snp, const size_t node_idx)
{
    TreeNode* const ref_node      = tree.node_ptr(node_idx);
    const TreeNode* comp_node     = tree.node_ptr(ref_node->root_idx);
    size_t          elements_used = 0;
   
    // Each of the snps which have been serached 
    for (size_t i = last_searched_snp; i >= 1; --i) {
        const size_t alignments       = tree.alignment_offsets[comp_node->haplo_idx + tree.reads];
        const size_t alignment_offset = tree.alignment_offsets[comp_node->haplo_idx];
        
        // Each of the alignments for a snp
        for (size_t j = alignment_offset; j < alignment_offset + alignments; ++j) {
            // If the aligned read crosses the snp index
            if (tree.read_info[tree.aligned_reads[j]].start_index() <= ref_node->haplo_idx &&
                tree.read_info[tree.aligned_reads[j]].end_index()   >= ref_node->haplo_idx ) {
                
                size_t row_offset  = tree.read_info[tree.aligned_reads[j]].offset();
                uint8_t ref_value  = tree.data[row_offset + ref_node->haplo_idx];
                uint8_t alignment  = tree.read_values[comp_node->align_idx + j];  
             
                if ((ref_value == ref_node->value && ref_value <= 1 && alignment == 1) ||
                    (ref_value != ref_node->value && ref_value <= 1 && alignment == 0)) {
                        // Optimal selection -- reduce the upper bound, dont increase lower bound
                        --ref_node->ubound; ++elements_used;
                } else if ((ref_value != ref_node->value && ref_value <= 1 && alignment == 1) ||
                           (ref_value == ref_node->value && ref_value <= 1 && alignment == 0)) {
                        // Non-optimal selection -- incread lower bound, don't reduce upper bound
                        ++ref_node->lbound; ++elements_used;
                }
            }
        }
        if (i > 1) {
            // Go back up the tree
            comp_node = tree.node_ptr(comp_node->root_idx);
        }
    }
    // For all remaining elements, we can reduce the upper bound
    ref_node->ubound -= (tree.snp_info[ref_node->haplo_idx].elements() - elements_used);
}

__global__ 
void map_level(internal::Tree tree          , BoundsGpu* snp_bounds      , const size_t last_searched_snp, 
               size_t* last_unaligned_idx   , size_t* alignment_offset   , const size_t prev_level_start , 
               const size_t this_level_start, const size_t nodes_in_level)
{
    // Set node parameters
    const size_t thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t node_idx  = thread_id + this_level_start;
    
    // If a valid thread
    if (thread_id < nodes_in_level) { 
        TreeNode* const         node       = tree.node_ptr(node_idx);
        const TreeNode* const   root_node  = tree.node_ptr(prev_level_start + (threadIdx.x / 2));
        const size_t last_unaligned_before = *last_unaligned_idx; 
       
        // Set some of the node parameters
        node->haplo_idx = snp_bounds[0].index;
        node->root_idx  = prev_level_start + (threadIdx.x / 2);
        node->lbound    = root_node->lbound;
        node->ubound    = root_node->ubound;
        node->value     = threadIdx.x % 2 == 0 ? 0 : 1;

        // Update the bounds for the node 
        map_leaf_bounds(tree, last_searched_snp, node_idx);

        if (thread_id == 0) {
            // Add the alignments to the overall alignments
            update_global_alignments(tree, last_unaligned_idx, node->haplo_idx);
            if (*last_unaligned_idx != last_unaligned_before) {
                // Set the start index in the alignment array
                tree.alignment_offsets[node->haplo_idx] = last_unaligned_before + 1; 
                tree.alignment_offsets[node->haplo_idx + tree.reads] = *last_unaligned_idx
                                                                     - last_unaligned_before;
            }
#ifdef DEBUG
            printf("ML : ALIGNEMENTS : %i\n", *last_unaligned_idx);
            for (size_t i = 0; i < tree.reads; ++i) 
                printf("%i ", tree.aligned_reads[i]);
            printf("\n");
#endif 
        }
        __syncthreads();
        
        // Add the alignments for the node 
        add_node_alignments(tree, last_unaligned_idx, node_idx, alignment_offset);
        if (thread_id == 0) {
            // Move the alignment offset for the next iteration
            *alignment_offset += (nodes_in_level * tree.alignment_offsets[node->haplo_idx + tree.reads]);
#ifdef DEBUG 
            printf("ML : ALIGN OFFSET : %i\n", *alignment_offset);
#endif
        }       
    }
}

__global__
void map_root_node(internal::Tree tree       , BoundsGpu* snp_bounds    , const size_t last_searched_snp,
                   size_t* last_unaligned_idx, size_t* alignment_offset , const size_t start_ubound     , 
                   const size_t device_index )
{
#ifdef DEBUG
    printf("Device Index : %i\n", device_index);
    printf("Start Bound  : %i\n" , start_ubound);

    for (size_t i = 0; i < tree.reads; ++i) printf("%i ", tree.aligned_reads[i]);
    printf("\n");
    for (size_t i = 0; i < tree.snps;  ++i) printf("%i ", tree.search_snps[i]);
    printf("\n");
#endif
    size_t last_unaligned_before = *last_unaligned_idx;
    
    TreeNode* const node  = tree.node_ptr(0);
    node->haplo_idx  = tree.search_snps[last_searched_snp];
    node->value      = 0; 
    node->lbound     = 0; 
    node->ubound     = start_ubound - tree.snp_info[node->haplo_idx].elements();
   
    // Add the alignments to the overall alignments
    update_global_alignments(tree, last_unaligned_idx, 0);
    if (*last_unaligned_idx != last_unaligned_before) {
        tree.alignment_offsets[0] = last_unaligned_before + 1; // Set the start index in the alignment array 
        tree.alignment_offsets[tree.reads] = *last_unaligned_idx - last_unaligned_before;
    }
    
    // Update the alignments values (how the reads are aligned) for this specific node
    add_node_alignments(tree, last_unaligned_idx, 0, alignment_offset);
    *alignment_offset += tree.alignment_offsets[tree.reads];        
    
#ifdef DEBUG
    printf("MRN : ALIGNMENTS : %i\n", tree.alignment_offsets[0]);
    printf("MRN : ALIGNMENTS : %i\n", *last_unaligned_idx);
    printf("MRN : ALIGN IDX  : %i\n", node->align_idx);
    for (size_t i = 0; i < tree.reads; ++i) 
        printf("%i ", tree.aligned_reads[i]);
    printf("\n");
    for (size_t i = 0; i < *last_unaligned_idx; ++i) 
        printf("%i ", tree.read_values[i]);
    printf("\n");
    printf("MRN : UPPER BOUND : %i\n", node->ubound);
#endif
    
}

}               // End namespace haplo
