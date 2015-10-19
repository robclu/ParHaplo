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

// Maps all the unsearched snps to a n array of BoundsGpu structs which can then be reduces
__global__ 
void map_unsearched_snps(internal::Tree tree, BoundsGpu* snp_bounds, const size_t last_searched_snp,
                         const size_t comp_node_idx)
{
    size_t read_start = 0, read_end = 0, row_offset = 0, opp = 0, same = 0;
    const size_t thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Get a reference to the bound parameter
    BoundsGpu* bounds = &snp_bounds[thread_id];
    
    // Index of the haplotype for the comparison
    const size_t ref_haplo_idx  = tree.search_snps[thread_id + last_searched_snp + 1];
    bounds->index               = ref_haplo_idx;
    bounds->offset              = thread_id + last_searched_snp + 1;
            
    // Get a pointer to the node
    const TreeNode* comp_node = &tree.nodes[comp_node_idx];
    
    // Go back up the tree and determine the 
    for (size_t i = last_searched_snp + 1; i > 0; --i) {
        // Go through the alignments for the node
        for (size_t j = 0; j < comp_node->alignments; ++j) {
            read_start = tree.read_info[comp_node->read_ids[j]].start_index();
            read_end   = tree.read_info[comp_node->read_ids[j]].end_index();
            // If any of the comp and ref snps (nodes) have overlapping positions
            if (read_start <= comp_node->haplo_idx && read_start <= ref_haplo_idx &&
                read_end   >= comp_node->haplo_idx && read_end   >= ref_haplo_idx ) {
                    row_offset = tree.read_info[comp_node->read_ids[j]].offset();

                    uint8_t comp_value = tree.data[row_offset + (comp_node->haplo_idx - read_start)];
                    uint8_t ref_value  = tree.data[row_offset + (ref_haplo_idx - read_start)];
                
                    if (comp_value == ref_value && ref_value <= 1) {
                        ++same;
                    } else if (comp_value != ref_value && ref_value <= 1 && comp_value <= 1) {
                        ++opp;
                    }
            }
        }
        if (i > 1) {
            // Move the pointer up the tree (towards the root)
            comp_node = tree.node_ptr(comp_node->root_idx);
        }
    }
    // Format the bounds 
    bounds->lower = min((unsigned int)same, (unsigned int)opp);
    bounds->upper = max((unsigned int)same, (unsigned int)opp);
    bounds->diff  = bounds->upper - bounds->lower;
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

// Very similar to the above function, but does the mapping for each left node of a level
__device__ 
void map_leaf_bounds(internal::Tree& tree, const size_t last_searched_snp, const size_t node_idx)
{
    size_t          elements_used = 0, row_offset = 0;
    uint8_t         ref_value     = 0;
    TreeNode* const ref_node      = tree.node_ptr(node_idx);
    const TreeNode* comp_node     = &tree.nodes[ref_node->root_idx];
    
    for (size_t i = last_searched_snp; i > 0; --i) {
        for (size_t j = 0; j < comp_node->alignments; ++j) {
            row_offset = tree.read_info[comp_node->read_ids[j]].offset();
            ref_value  = tree.data[row_offset + ref_node->haplo_idx];
             
            if ((ref_value == ref_node->value && ref_value <= 1 && comp_node->read_ids[j] == 1) ||
                (ref_value != ref_node->value && ref_value <= 1 && comp_node->read_ids[j] == 0)) {
                    // Optimal selection -- reduce the upper bound, dont increase lower bound
                    --ref_node->ubound; ++elements_used;
            } else if ((ref_value != ref_node->value && ref_value <= 1 && comp_node->read_ids[j] == 1) ||
                       (ref_value == ref_node->value && ref_value <= 1 && comp_node->read_ids[j] == 0)) {
                    // Non-optimal selection -- incread lower bound, don't reduce upper bound
                    ++ref_node->lbound; ++elements_used;
            }
        }
        if (i > 1) {
            comp_node = tree.node_ptr(comp_node->root_idx);
        }
    }
    // For all remaining elements, we can reduce the upper bound
    ref_node->ubound -= (tree.snp_info[ref_node->haplo_idx].elements() - elements_used);
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
    const size_t reductions     = (size_t)(ceil(log2((double)elements)));
    size_t reduction_threads    = elements - 1;
    size_t snp_idx_other        = 0;
    
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
        // If we came from an odd number of bounds, the last one just fetches a value
        if (reduction_threads % 2 == 1) {
            if (thread_id == (reduction_threads / 2)) {
                // There is an odd number of elements in the array,
                // The last thread just needs to move a value 
                BoundsGpu temp          = snp_bounds[thread_id];
                snp_bounds[thread_id]  = snp_bounds[snp_idx_other];
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


// Add any unaligned reads to a node (snp site)
__device__ 
void add_alignments(internal::Tree& tree, size_t* last_unaligned_idx, 
                    const size_t node_idx, const size_t thread_id)
{
    size_t read_offset = 0, align_count = 0, element_position = 0;
    
    // Get a pointer to the node
    TreeNode* node = tree.node_ptr(node_idx);
    
    // Allocate some temp space for the node alignments (this is usually small)
    size_t* alignments = new size_t[tree.snp_info[node->haplo_idx].elements()]; 
    size_t* values     = new size_t[tree.snp_info[node->haplo_idx].elements()];
    size_t* indices;
    
    if (thread_id == 0) {
        indices = new size_t[tree.snp_info[node->haplo_idx].elements()];
    }          
    
    // Go through all the unaligned reads
    for (size_t i = *last_unaligned_idx; i < tree.reads; ++i) {
        if (tree.aligned_reads[i] >= tree.snp_info[node->haplo_idx].start_index() &&
            tree.aligned_reads[i] <= tree.snp_info[node->haplo_idx].end_index()  ) {
        
            // If the row crosses the snp
            if (tree.read_info[tree.aligned_reads[i]].start_index() <= node->haplo_idx &&
                tree.read_info[tree.aligned_reads[i]].end_index()   >= node->haplo_idx ) {
            
                // Get the offset in memory of the start of the read
                read_offset         = tree.read_info[tree.aligned_reads[i]].offset();
                element_position    = node->haplo_idx - tree.read_info[tree.aligned_reads[i]].start_index();
                uint8_t  element    = tree.data[read_offset + element_position];
            
                // Do the alignement for the reads to this node
                if ((element == 0 && node->value == 0) || (element == 1 && node->value == 1)) {
                    values[align_count] = 1; alignments[align_count++] = tree.aligned_reads[i];
                    if (thread_id == 0) indices[align_count - 1] = i; 
                } else if ((element == 0 && node->value == 1) || (element == 1 && node->value == 0)) {
                    values[align_count] = 0; alignments[align_count++] = tree.aligned_reads[i];
                    if (thread_id == 0) indices[align_count - 1] = i; 
                }
            }
        }
    }
   
    // Move the found alignments to the node
    node->alignments     = align_count;
    node->read_ids       = new size_t[align_count];
    node->read_values    = new uint8_t[align_count];
   
    if (thread_id == 0) node->indices = new size_t[align_count];
    
    for (size_t i = 0; i < align_count; ++i) {
        node->read_ids[i]    = alignments[i];
        node->read_values[i] = (uint8_t)values[i];
        if (thread_id == 0) node->indices[i] = indices[i];
    }
    // Clean memory
    free(alignments); free(values);
    if (thread_id == 0) free(indices);
}

// Updates the array of reads which have been aligned
__device__
void update_global_alignments(internal::Tree& tree, size_t* last_unaligned_idx, const size_t node_idx) 
{
    const TreeNode* const node = tree.node_ptr(node_idx);
    
    for (size_t i = *last_unaligned_idx; i < *last_unaligned_idx + node->alignments; ++i) {
        size_t temp = tree.aligned_reads[node->indices[i - *last_unaligned_idx]];
        tree.aligned_reads[node->indices[i - *last_unaligned_idx]] = tree.aligned_reads[i];
        tree.aligned_reads[i] = temp;
    } 
    *last_unaligned_idx += node->alignments;
}

__global__ 
void map_level(internal::Tree tree, BoundsGpu* snp_bounds, const size_t last_searched_snp, 
               size_t* last_unaligned_idx, const size_t prev_level_start, const size_t   this_level_start)
{
    // Set node parameters
    const size_t            node_idx  = threadIdx.x + this_level_start;
    TreeNode* const         node      = tree.node_ptr(node_idx);
    const TreeNode* const   root_node = tree.node_ptr(prev_level_start + (threadIdx.x / 2));
   
    // Set some of the node parameters
    node->haplo_idx = snp_bounds[0].index;
    node->root_idx  = prev_level_start + (threadIdx.x / 2);
    node->lbound    = root_node->lbound;
    node->ubound    = root_node->ubound;
    node->value     = threadIdx.x % 2 == 0 ? 0 : 1;
    
    // Update the bounds for the node 
    map_leaf_bounds(tree, last_searched_snp, node_idx);

    // Add the alignments for the node 
    add_alignments(tree, last_unaligned_idx, node_idx, threadIdx.x);
    
    // The first thread needs to update the global alignments as well
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        update_global_alignments(tree, last_unaligned_idx, this_level_start);
        
#ifdef DEBUG
        printf("ML : ALIGNEMENTS : %i\n", *last_unaligned_idx);
        for (size_t i = 0; i < tree.reads; ++i) 
            printf("%i ", tree.aligned_reads[i]);
        printf("\n");
#endif
    }   
    __syncthreads();
}

// Swaps two nodes in the node array
__device__ 
void swap(TreeNode* left, TreeNode* right) 
{
    TreeNode temp = *left;
    *left         = *right;
    *right        = temp;
}

// Compares two tree nodes 

// "Reduces" the nodes, removing all the nodes with a lower bounds > than the highest upper bound
// by movind them to the end of the array and then returning the index of the end of the array
// also sorts the nodes by lower bound -- requires 2*log_2(N) operations
__global__
void reduce_level(internal::Tree tree, size_t* start_node_idx, const size_t* num_nodes)
{
    const size_t thread_id          = blockIdx.x * blockDim.x + threadIdx.x;
    size_t       reduction_threads  = *num_nodes - 1;
    size_t       node_idx           = 0;
    
    while (reduction_threads > 1) {
        // ----------------------------------- SORT STEP ----------------------------------------------------
        node_idx = thread_id + *start_node_idx;
        if (thread_id % 2 == 0 && thread_id != reduction_threads) {
            if (tree.node_ptr(node_idx)->lbound > tree.node_ptr(node_idx + 1)->lbound) {
                swap(tree.node_ptr(node_idx), tree.node_ptr(node_idx + 1));
            } 
            // node_idx now has the lower bound 
            if (tree.node_ptr(node_idx)->ubound <= tree.node_ptr(node_idx + 1)->lbound) {
                // Add to counters here 
            }
        } 
        __synthreads();
        
        // ----------------------------------- SWAP STEP ----------------------------------------------------
        
        if (reduction_threads % 2 == 1 && thread_id != reduction_threads) {
            if (tree.node_ptr(node_idx)->lbound > tree.node_ptr(node_idx + 1)->lbound) {
                swap(tree.node_ptr(node_idx), tree.node_ptr(node_idx + 1));
            } 
            // node_idx is lower bound again
            if (tree.node_ptr(node_idx)->ubound <= tree.node_ptr(node_idx + 1)->lbound) {
                // Add to the counters again
            }
        }
        // Modify number of threads for next iteration
        reduction_threads /= 2;
        if (reduction_threads % 2 == 1) reduction_threads += 1;
        __syncthreads();
    }
    
}

__global__
void map_root_node(internal::Tree tree, BoundsGpu* snp_bounds, const size_t last_searched_snp,
                   size_t* last_unaligned_idx, const size_t start_ubound, const size_t device_index)
{
    // DEBUG 
    printf("Device Index : %i\n", device_index);
    printf("Start Bound  : %i\n" , start_ubound);
    
    struct cudaDeviceProp device_properties;                // So that we can know the max number of threads

    // ---------------------------------------- ROOT NODE -------------------------------------------------

    // DEBUG
    for (size_t i = 0; i < tree.reads; ++i) printf("%i\n", tree.aligned_reads[i]);
    for (size_t i = 0; i < tree.snps;  ++i) printf("%i\n", tree.search_snps[i]);
    
    TreeNode& node = tree.nodes[0];
    node.haplo_idx = tree.search_snps[last_searched_snp];
    node.node_idx  = 0; node.value  = 0;
    
    // Set the alignments for the tree root
    add_alignments(tree, last_unaligned_idx, 0, 0);
   
    // Add the alignments to the overall alignments
    update_global_alignments(tree, last_unaligned_idx, 0);
   
    // Update the bounds for the tree
    node.lbound = 0; node.ubound = start_ubound - tree.snp_info[node.haplo_idx].elements();
    
#ifdef DEBUG
    printf("MRN : ALIGNMENTS : %i\n", *last_unaligned_idx);
    for (size_t i = 0; i < tree.reads; ++i) 
        printf("%i ", tree.aligned_reads[i]);
    printf("\n");
    printf("MRN : UPPER BOUND : %i\n", node.ubound);
#endif
    
}

}               // End namespace haplo
