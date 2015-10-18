#include "cuda.h"
#include "math.h"
#include "tree_internal.h"

namespace haplo {

struct BoundsGpu {

    size_t lower;               // The lower value for the bound    
    size_t upper;               // The upper value for the bound
    size_t diff;                // The difference between the upper and lower bound
    size_t index;               // The snp index the bound represents
    
    // Overlaod the equality operator
    __device__
    void operator=(BoundsGpu& other) 
    {
        lower           = other.lower;
        upper           = other.upper;
        diff            = other.diff;
        index           = other.index;
    }
};

// Maps all the unsearched snps to a n array of BoundsGpu structs which can then be reduces
__global__ 
void map_unsearched_snps(internal::Tree tree, BoundsGpu* snp_bounds, const size_t comp_node_idx)
{
    size_t read_start = 0, read_end = 0, row_offset = 0, opp = 0, same = 0;
    
    const size_t thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Get a reference to the bound parameter
    BoundsGpu* bounds = &snp_bounds[thread_id];
    
    // Index of the haplotype for the comparison
    const size_t ref_haplo_idx = tree.search_snps[thread_id + tree.last_searched_snp + 1];
    bounds->index              = ref_haplo_idx;
            
    // Get a pointer to the node
    const TreeNode* comp_node = &tree.nodes[comp_node_idx];
    
    // Go back up the tree and determine the 
    for (size_t i = tree.last_searched_snp + 1; i > 0; --i) {
         
        // DEBUG  
        printf("Node HI : %i ,", comp_node->haplo_idx);
        printf("Ref HI : %i\n" , ref_haplo_idx);
        
        // Go through the alignments for the node
        for (size_t j = 0; j < comp_node->alignments; ++j) {
            read_start = tree.read_info[comp_node->read_ids[j]].start_index();
            read_end   = tree.read_info[comp_node->read_ids[j]].end_index();
            
            // If any of the comp and ref snps (nodes) have overlapping positions
            if (read_start <= comp_node->haplo_idx && read_start <= ref_haplo_idx &&
                read_end   >= comp_node->haplo_idx && read_end   >= ref_haplo_idx ) {
                row_offset = tree.read_info[comp_node->read_ids[j]].offset();

                uint8_t comp_value = tree.data[row_offset + comp_node->haplo_idx];
                uint8_t ref_value  = tree.data[row_offset + ref_haplo_idx];
                
                if (comp_value == ref_value && ref_value <= 1) {
                    ++same;
                    // DEBUG
        //            printf("Adding Same : %i, %i\n", comp_node->read_ids[j], ref_haplo_idx);
                } else if (comp_value != ref_value && ref_value <= 1 && comp_value <=1) {
                    ++opp;
                    // DEBUG
        //            printf("Adding Oppp : %i, %i\n", comp_node->read_ids[j], ref_haplo_idx);
                }
            }
        }
        if (i > 1) {
            // Move the pointer up the tree (towards the root)
            const size_t root_idx = comp_node->root_idx;
            comp_node = tree.node_ptr(root_idx);
            
            // DEBUG
          //  printf("Reassign : Root ID : %i : Node ID :\n", root_idx, comp_node->node_idx); 
        }
    }
    // Format the bounds 
    bounds->lower = min((unsigned int)same, (unsigned int)opp);
    bounds->upper = max((unsigned int)same, (unsigned int)opp);
    bounds->diff  = bounds->upper - bounds->lower;
    
    // DEBUG
    //printf("Bounds Diff Inline: %i\n", bounds->diff);
}

// Very similar to the above function, but does the mapping for each left node of a level
__device__ 
void map_leaf_bounds(internal::Tree& tree, const size_t node_idx)
{
    size_t          elements_used = 0, row_offset = 0;
    uint8_t         ref_value     = 0;
    TreeNode* const ref_node      = tree.node_ptr(node_idx);
    const TreeNode* comp_node     = &tree.nodes[ref_node->root_idx];
    
    for (size_t i = tree.last_searched_snp; i > 0; --i) {
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
            
            // DEBUG 
            printf("Moving Up Tree\n");
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

// "Reduce" function for the list of selection params to determine the best param 
//  The reduction is done such that the resulting bound has the highest diff and 
//  lowest index (this is to mean that the snp is most correlated to the snps already
//  searched and is also closest to them (hence the lowest index))
__global__ 
void reduce_unsearched_snps(BoundsGpu* snp_bounds, const size_t elements)
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
}


// Add any unaligned reads to a node (snp site)
__device__ 
void add_alignments(internal::Tree& tree, const size_t node_idx, const size_t thread_id)
{
    size_t read_offset = 0, align_count = 0, element_position = 0;
    
    // Get a pointer to the node
    TreeNode* node = tree.node_ptr(node_idx);
   
    // DEBUG
    //printf("AlignS : %i\n", node->haplo_idx);
    
    // Allocate some temp space for the node alignments (this is usually small)
    size_t* alignments = new size_t[tree.snp_info[node->haplo_idx].elements()]; 
    size_t* values     = new size_t[tree.snp_info[node->haplo_idx].elements()];
    size_t* indices;
    
    if (thread_id == 0) {
        indices = new size_t[tree.snp_info[node->haplo_idx].elements()];
    }
                
    // Go through all the unaligned reads
    for (size_t i = tree.last_unaligned_idx; i < tree.reads; ++i) {
        if (tree.aligned_reads[i] >= tree.snp_info[node->haplo_idx].start_index() &&
            tree.aligned_reads[i] <= tree.snp_info[node->haplo_idx].end_index()  ) {
        
            // If the row crosses the snp
            if (tree.read_info[tree.aligned_reads[i]].start_index() <= node->haplo_idx &&
                tree.read_info[tree.aligned_reads[i]].end_index()   >= node->haplo_idx ) {
            
                // Get the offset in memory of the start of the read
                read_offset         = tree.read_info[tree.aligned_reads[i]].offset();
                element_position    = node->haplo_idx - tree.read_info[tree.aligned_reads[i]].start_index();
                auto element = tree.data[read_offset + element_position];
            
                // Do the alignement for the reads to this node
                if ((element == 0 && node->value == 0) || (element == 1 && node->value == 1)) {
                    // DEBUG
     //               printf("Align1 : %i\n", i);
                    values[align_count] = 1; alignments[align_count++] = tree.aligned_reads[i];
                    if (thread_id == 0) indices[align_count - 1] = i; 
                } else if ((element == 0 && node->value == 1) || (element == 1 && node->value == 0)) {
                    values[align_count] = 0; alignments[align_count++] = tree.aligned_reads[i];
                    if (thread_id == 0) indices[align_count - 1] = i; 
                    // DEBUG
     //               printf("Align2 : %i\n", i);
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


__global__ 
void map_level(internal::Tree tree            , BoundsGpu* snp_bounds, const size_t prev_level_start, 
               const size_t   this_level_start)
{
    // Set node parameters
    size_t      node_idx  = threadIdx.x + this_level_start;
    TreeNode*   node      = tree.node_ptr(node_idx);
    TreeNode*   root_node = tree.node_ptr(prev_level_start + (threadIdx.x / 2));
   
    printf("Reached!\n");
    
    // Set some of the node parameters
    node->haplo_idx = snp_bounds[0].index;
    node->root_idx  = prev_level_start + (threadIdx.x / 2);
    node->lbound    = root_node->lbound;
    node->ubound    = root_node->ubound;
    node->value     = threadIdx.x % 2 == 0 ? 0 : 1;
  
    // DEBUG
   // printf("H ID : %i\n", node->haplo_idx);
    
    // DEBUG
    printf("Lbound : %i " , node->lbound); 
    printf("Ubound : %i\n", node->ubound); 
    
    // Update the bounds for the node 
    map_leaf_bounds(tree, node_idx);
  
    // DEBUG
    printf("Lbound : %i " , node->lbound); 
    printf("Ubound : %i\n", node->ubound); 

    // Add the alignments for the node 
    add_alignments(tree, node_idx, threadIdx.x);
    
    // DEBUG
    printf("Alignments : %i\n", node->alignments);
}

__device__ 
void swap_search_snp_indices(internal::Tree& tree, const size_t swap_idx)
{
    ++tree.last_searched_snp;
    if (swap_idx != tree.last_searched_snp) {
        const size_t temp                        = tree.search_snps[tree.last_searched_snp];
        tree.search_snps[tree.last_searched_snp] = tree.search_snps[swap_idx];
        tree.search_snps[swap_idx]               = temp;
    }
}


// Updates the array of reads which have been aligned
__device__
void update_global_alignments(internal::Tree& tree, const size_t node_idx) 
{
    // DEBUG 
    printf("Node ID Align : %i\n", node_idx);
    
    const TreeNode* const node = tree.node_ptr(node_idx);
   
    printf("Aligns: %i\n", node->alignments);
    
    for (size_t i = tree.last_unaligned_idx; i < tree.last_unaligned_idx + node->alignments; ++i) {
        size_t temp = tree.aligned_reads[node->indices[i - tree.last_unaligned_idx]];
        tree.aligned_reads[node->indices[i - tree.last_unaligned_idx]] = tree.aligned_reads[i];
        tree.aligned_reads[i] = temp;
    } 
    tree.last_unaligned_idx += node->alignments;
    printf("Finished Align!\n");
}


__device__ size_t       prev_level_start = 0;           // The index of the first nod in the previous level
__device__ size_t       this_level_start = 1;           // For each level, this is the index in the 
                                                        // node array of the first element in the level 
__device__ size_t       nodes_in_level   = 2;           // The number of nodes (sub-branches) in the level
__device__ size_t       unsearched_snps;                // The number of unsearched snps
__device__ size_t       comp_node_idx;                  // The index of the comparison node for node seletion

__global__ void search_tree(internal::Tree tree, BoundsGpu* snp_bounds, size_t start_ubound, size_t device_index)
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
    node.haplo_idx = tree.search_snps[tree.last_searched_snp];
    node.node_idx  = 0; node.value  = 0;
    
    // Set the alignments for the tree root
    add_alignments(tree, 0, 0);

    // Add the alignments to the overall alignments
    update_global_alignments(tree, 0);
    
    // DEBUG
    for (size_t i = 0; i < tree.last_unaligned_idx; ++i) {
        printf("%i ", tree.aligned_reads[i]);
    } printf("\n");

    // Go over all the nodes that have not been searched and see how correlated they are
    comp_node_idx   = node.node_idx; 
    unsearched_snps = tree.snps - tree.last_searched_snp - 1;

    // Perform a "mapping" step to map all the unsearched snps to their potenetial bounds
    map_unsearched_snps<<<1, unsearched_snps>>>(tree, snp_bounds, comp_node_idx);
    if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for SNP Map!\n");
    cudaDeviceSynchronize(); 
    __syncthreads();
    
    // Do the reduction to get the index of the next node to search
    reduce_unsearched_snps<<<1, unsearched_snps>>>(snp_bounds, unsearched_snps);
    cudaDeviceSynchronize();
    __syncthreads();
    
    --unsearched_snps; // We have searched another snp
    
    // DEBUG 
    for (size_t i = 0; i < unsearched_snps; ++i) {
        printf("%i ", snp_bounds[i].diff);
    } printf("\n");
    for (size_t i = 0; i < unsearched_snps; ++i) {
        printf("%i ", snp_bounds[i].index);
    } printf("\n");
    
    // Reduction has found the next node, modify the searned snp array
    swap_search_snp_indices(tree, snp_bounds[0].index);
    
    // Update the bounds for the tree
    node.lbound = 0; node.ubound = start_ubound - tree.snp_info[node.haplo_idx].elements();
    
    // ----------------------------------- STARTING TREE SEARCH -------------------------------------------------
   
    // Iterate through level
    size_t iteration = 0;
    while (tree.last_searched_snp < tree.snps && counter++ < 13) {
        if (cudaSuccess != cudaGetLastError()) printf("Error!\n");
       
        printf("Starting new iteration!\n");
        
        // We can only launch 
        // Perform a "mapping" step, which maps the bounds onto the nodes in the level
        map_level<<<1, nodes_in_level>>>(tree, snp_bounds, prev_level_start, this_level_start);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for Level Map!\n");
        cudaDeviceSynchronize();
        __syncthreads();
        
        // Add the new alignments to the array of overall aligned reads -- use first node
        if (threadIdx.x == 0) {
            update_global_alignments(tree, this_level_start);
        }
        __syncthreads(); 
        
        // Now "reduce" the level, which is essentially a pruning step
        
        // Map unsearched snps
        map_unsearched_snps<<<1, unsearched_snps>>>(tree, snp_bounds, this_level_start);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for SNP Map!\n");
        cudaDeviceSynchronize();
        __syncthreads();
       
        // DEBUG
        printf("Before\n");
        for (size_t i = 0; i < unsearched_snps; ++i) {
            printf("%i ", snp_bounds[i].diff);
        } printf("\n");
        for (size_t i = 0; i < unsearched_snps; ++i) {
            printf("%i ", snp_bounds[i].index);
        } printf("\n");
            
        reduce_unsearched_snps<<<1, unsearched_snps>>>(snp_bounds, unsearched_snps);
        if (cudaSuccess != cudaGetLastError()) printf("Kernel Launch Error for SNP Reduce!\n");
        cudaDeviceSynchronize();
     
        if (threadIdx.x == 0) {   
            --unsearched_snps;
        }
        __syncthreads();
        
        // DEBUG
        for (size_t i = 0; i < unsearched_snps; ++i) {
            printf("%i ", snp_bounds[i].diff);
        } printf("\n");
        for (size_t i = 0; i < unsearched_snps; ++i) {
            printf("%i ", snp_bounds[i].index);
        } printf("\n");

        // DEBUG
        printf("Aligned Reads : %i\n", tree.last_unaligned_idx);
        for (size_t i = 0; i < tree.reads; ++i) {
            printf("%i ", tree.aligned_reads[i]);
        } printf("\n");
       
        printf("Lowest Index :%i\n", snp_bounds[0].index);
        
        // Swap the indices of the last unsearched and the about to be searched snps
        swap_search_snp_indices(tree, snp_bounds[0].index);
        
        // DEBUG
        printf("Searched SNPS :\n");
        for (size_t i = 0; i < tree.snps; ++i) {
            printf("%i ", tree.search_snps[i]);
        } printf("\n");
        
        // Change the indices of the start and end of the level
        prev_level_start  = this_level_start;
        this_level_start += nodes_in_level;
        nodes_in_level   *= 2;
        
        printf("Nodes in level : %i\n", nodes_in_level);
        printf("Start Node     : %i\n", this_level_start);
        printf("Last Node      : %i\n", this_level_start + nodes_in_level);
        printf("Last SNP       : %i\n", tree.last_searched_snp);
    }
        
        
    /*
    size_t max = 0, index = 0; BoundsGpu bounds_temp, bounds_final;
    for (size_t i = tree.last_searched_snp + 1; i < tree.snps; ++i) {yy
        bounds_temp = compare_snps(tree, node.node_idx, i);
        if (bounds_temp.diff > max) { 
            index        = i; 
            max          = bounds_temp.diff; 
            bounds_final = bounds_temp;
        }
        printf("Result : %i\n", bounds_temp.diff);
        printf("Max    : %i\n", max   );
        printf("Index  : %i\n", index );
   }
  
    // DEBUG
    printf("Most Correlated : %i ", index);
    printf("Correlation : %i\n"   , max  );
    
    // The first node has now been searched
    tree.last_searched_snp++;
    
    // Swap the SNPs here
   
    // Set the upper and lower bounds based on the most correlated node
    node.lbound += bounds_final.lower; node.ubound -= bounds_final.upper;
    
   */
/*    
    // Make the next 2 nodes point back to this one
    TreeNode& left_child = tree.nodes[1]; TreeNode& right_child = tree.nodes[2];
    left_child.root_idx  = 0; right_child.root_idx  = 0;  
    left_child.value     = 1; right_child.value     = 1;
    left_child.node_idx  = 1; right_child.node_idx  = 2;
        
    printf("%i rest\n ",*result);
    
    left_child.haplo_idx = 1; right_child.haplo_idx  = 1;
    printf("%i rest\n ", left_child.node_idx);
    
    add_alignments(tree, tree.read_info, tree.reads, tree.snp_info, tree.data, &left_child , aligned, tree.last_aligned);
    add_alignments(tree, tree.read_info, tree.reads, tree.snp_info, tree.data, &right_child, aligned, tree.last_aligned);
   
    printf("%i rest\n ", left_child.node_idx);

    //left_child.lbound    = bounds.same; left_child.ubound  = *ubound - (bounds.same + left_child.alignments);
    //right_child.lbound   = bounds.opps; right_child.ubound = *ubound - (bounds.opps + right_child.alignments);
    
    // Add the alignments due to these nodes -- same if we use right or left child
    for (size_t i = tree.last_aligned; i < tree.last_aligned + left_child.alignments; ++i) {
            aligned[left_child.read_ids[i - tree.last_aligned]] = aligned[i];
            aligned[i] = left_child.read_ids[i - tree.last_aligned];
    } tree.last_aligned += left_child.alignments;

    printf("%i rest\n ", left_child.node_idx);

    temp = 0;
    for (size_t i = tree.last_searched_snp + 1; i < tree.snps; ++i) {
        res = compare_snps(tree, tree.read_info, tree.data, &left_child, i, &tree.last_searched_snp);
        if (res > temp) { *result = i; temp = res; }
    }

    printf("%i : %i rest\n ", temp, *result);

//    for (size_t i = 0; i < tree.snps; ++i) 
//        printf("%i ", tree.search_snps[i]);

    printf("\nres: %i, %i\n", right_child.lbound, left_child.node_idx);
    
    //start_index = 1; snp_end = tree.last_searched_snp; haplo_idx = *result; snps = tree.snps; reads = tree.reads;
    //search_helper<<<2, 1>>>(tree, tree.node_manager.nodes, tree.read_info, tree.snp_info, tree.data, &haplo_idx, &snp_end, &start_index, &snps, &reads);
 
    // Now we can start the iterative search 
    
    //size_t a = compare_snps(tree, &left_child, 2);
    tree.last_searched_snp++;
   
   // *result = left_child.alignments;
    
    //free(aligned);
    */
}

}               // End namespace haplo
