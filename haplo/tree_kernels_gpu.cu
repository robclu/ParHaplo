#include "cuda.h"
#include "math.h"
#include "tree_internal.h"

namespace haplo {

struct BoundsGpu {

bool   lower_is_same;
size_t lower;
size_t upper;
size_t diff;

};
__device__ 
BoundsGpu compare_snps(internal::Tree& tree, const size_t comp_node_idx, const size_t ref_haplo_idx)
{
//    size_t node_idx = node.node_idx;
    size_t read_start = 0, read_end     = 0, row_offset  = 0,
           opp        = 0, same         = 0;
    
    // The bounds for this instance
    BoundsGpu bounds;
            
    // Get a pointer to the node
    const TreeNode* comp_node = &tree.nodes[comp_node_idx];
    
    // Go back up the tree and determine the 
    for (size_t i = tree.last_searched_snp + 1; i > 0; --i) {
      
        // DEBUG  
        printf("Node ID : %i Node HI : %i, Ref HI : %i\n", comp_node->node_idx, comp_node->haplo_idx, ref_haplo_idx);
        
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
                    printf("Adding Same : %i, %i\n", comp_node->read_ids[j], ref_haplo_idx);
                } else if (comp_value != ref_value && ref_value <= 1 && comp_value <=1) {
                    ++opp;
                    // DEBUG
                    printf("Adding Oppp : %i, %i\n", comp_node->read_ids[j], ref_haplo_idx);
                }
            }
        }
        if (i > 1) {
            // Move the pointer up the tree (towards the root)
            const size_t root_idx = comp_node->root_idx;
            comp_node = tree.node_ptr(root_idx);
            
            // DEBUG
            printf("Reassign : Root ID : %i : Node ID :\n", root_idx, comp_node->node_idx); 
        }
    }
    // Format the bounds 
    bounds.lower = min((unsigned int)same, (unsigned int)opp);
    bounds.upper = max((unsigned int)same, (unsigned int)opp);
    bounds.diff  = bounds.upper - bounds.lower;
    bounds.lower_is_same = same <= opp ? true : false;
    
    return bounds;
}

__device__ size_t* realloc(size_t& old_size, size_t new_size, size_t* old)
{
    size_t* new_array = (size_t*)malloc(new_size * sizeof(size_t));

    for (size_t i = 0; i < old_size; i++) new_array[i] = old[i];

    free(old);
    return new_array;
}

// Add any unaligned reads to a node (snp site)
__device__ 
void add_alignments(internal::Tree& tree, const size_t node_idx)
{
    size_t read_offset = 0, align_count = 0;
    
    // Get a pointer to the node
    TreeNode* node = tree.node_ptr(node_idx);
   
    // DEBUG
    printf("AlignS : %i\n", node->haplo_idx);
    
    // Allocate some temp space for the node alignments (this is usually small)
    size_t* alignments = new size_t[tree.snp_info[node->haplo_idx].end_index() - 
                                    tree.snp_info[node->haplo_idx].start_index()];
    size_t* values     = new size_t[tree.snp_info[node->haplo_idx].end_index() - 
                                    tree.snp_info[node->haplo_idx].start_index()];     
    
    // Go through all the unaligned reads
    for (size_t i = tree.last_unaligned_idx; i < tree.reads; ++i) {
        if (i >= tree.snp_info[node->haplo_idx].start_index() &&
            i <= tree.snp_info[node->haplo_idx].end_index()  ) {
           
            // DEBUG
            printf("Test\n");
            
            // If the row crosses the snp
            if (tree.read_info[tree.aligned_reads[i]].start_index() <= node->haplo_idx &&
                tree.read_info[tree.aligned_reads[i]].end_index()   >= node->haplo_idx ) {
            
                // Get the offset in memory of the start of the read
                read_offset  = tree.read_info[i].offset();
                auto element = tree.data[read_offset + node->haplo_idx];
            
                // Do the alignement for the reads to this node
                if ((element == 0 && node->value == 0) || (element == 1 && node->value == 1)) {
                    // DEBUG
                    printf("Align1 : %i\n", i);
                    values[align_count] = 1; alignments[align_count++] = tree.aligned_reads[i];
                } else if ((element == 0 && node->value == 1) || (element == 1 && node->value == 0)) {
                    values[align_count] = 0; alignments[align_count++] = tree.aligned_reads[i];
                    // DEBUG
                    printf("Align2 : %i\n", i);
                }
            }
        }
    }
   
    // Move the found alignments to the node
    node->alignments     = align_count;
    node->read_ids       = new size_t[align_count];
    node->read_values    = new uint8_t[align_count];
    
    for (size_t i = 0; i < align_count; ++i) {
        node->read_ids[i]    = alignments[i];
        node->read_values[i] = (uint8_t)values[i];
    }
    // Clean memory
    free(alignments); free(values);
}


__global__ void search_helper(internal::Tree tree, TreeNode* nodes, ReadInfo* read_info, SnpInfoGpu* snp_info, uint8_t* data, 
        size_t* haplo_idx, size_t* last_snp, size_t* start_idx, size_t* snps, size_t* reads)
{
    // Set node parameters
    BoundsGpu bounds;
    size_t      node_idx  = threadIdx.x + *start_idx;
    TreeNode*   node      = &nodes[node_idx];
    node->haplo_idx         = *haplo_idx;
    
    //size_t res = compare_snps(read_info, data, &node, 1, last_snp, &bounds);
   // size_t temp = 0, res = 0;
  //  for (size_t i = *last_snp + 1; i < *snps; ++i) {
  //      res = compare_snps(tree, read_info, data, node, 2, *last_snp, &bounds);
  //      if (res > temp) { node->haplo_idx = i; temp = res; }
   // }
    
    //add_alignments(tree, read_info, *reads, snp_info, tree.data, node , set_alignments, last_aligned);
}

__device__ size_t start_node_index = 1;                     // For each level, this is the index in the 
                                                            // node array of the first element in the level 
__device__ size_t nodes_in_level = 2;                       // The number of nodes (sub-branches) in the level

__global__ void search_tree(internal::Tree tree, size_t start_ubound, size_t device_index)
{
    // DEBUG 
    printf("Device Index : %i\n", device_index);
    printf("Start Bound  : %i\n" , start_ubound);
    
    struct cudaDeviceProp device_properties;                // So that we can know the max number of threads

    // Get the properties of the device 
    //
    //cudaError_t status = 
    // ---------------------------------------- ROOT NODE -------------------------------------------------

    // DEBUG
    for (size_t i = 0; i < tree.reads; ++i) printf("%i\n", tree.aligned_reads[i]);
    for (size_t i = 0; i < tree.snps;  ++i) printf("%i\n", tree.search_snps[i]);
    
    TreeNode& node = tree.nodes[0];
    node.haplo_idx = tree.search_snps[tree.last_searched_snp];
    node.node_idx  = 0; node.value  = 0;
    
    // Set the alignments for the tree root
    add_alignments(tree, 0);
    
    // DEBUG 
    printf("AlignF : %i\n", node.alignments);

    // Add the alignments to the overall alignments
    for (size_t i = tree.last_unaligned_idx; i < tree.last_unaligned_idx + node.alignments; ++i) {
        tree.aligned_reads[node.read_ids[i - tree.last_unaligned_idx]] = tree.aligned_reads[i];
        tree.aligned_reads[i] = node.read_ids[i - tree.last_unaligned_idx];
    } tree.last_unaligned_idx += node.alignments;
    
    // DEBUG 
    printf("AlignL : %i\n", tree.last_unaligned_idx);
    for (size_t i = 0; i < tree.last_unaligned_idx; ++i) printf("%i ", tree.aligned_reads[i]);
    printf("\nLast Searched: %i\n", tree.last_searched_snp);    

    // Go over all the nodes that have not been searched and see how correlated they are
    size_t max = 0, index = 0; BoundsGpu bounds_temp, bounds_final;
    for (size_t i = tree.last_searched_snp + 1; i < tree.snps; ++i) {
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
    printf("Most Correlated : %i Correlation : %i\n", index, max);
    
    // The first node has now been searched
    tree.last_searched_snp++;
    
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
