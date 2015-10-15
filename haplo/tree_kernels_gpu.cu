#include "cuda.h"
#include "math.h"
#include "tree_internal.h"

namespace haplo {

struct Bounds {

size_t same;
size_t opps;

};

__device__ size_t compare_snps(internal::Tree& tree, TreeNode* comp_node, const size_t ref_haplo_idx, 
        Bounds* bounds)
{
//    size_t node_idx = node.node_idx;
    size_t read_start = 0, read_end     = 0, row_offset  = 0,
           opp        = 0, same         = 0;

    // Go back up the tree and determine the 
    for (size_t i = tree.last_searched_snp + 1; i > 0; --i) {
        // Go through the alignments
        for (size_t j = 0; j < comp_node->alignments; ++j) {
            read_start = tree.read_info[comp_node->read_ids[j]].start_index();
            read_end   = tree.read_info[comp_node->read_ids[j]].end_index();
            
            if (read_start <= comp_node->haplo_idx && read_start <= ref_haplo_idx &&
                read_end   >= comp_node->haplo_idx && read_end   >= ref_haplo_idx ) {
                row_offset = tree.read_info[comp_node->read_ids[j]].offset();

                auto comp_value = tree.data[row_offset + comp_node->haplo_idx];
                auto ref_value  = tree.data[row_offset + ref_haplo_idx];
                
                if (comp_value == ref_value && ref_value <= 1) {
                    ++same;
                    comp_node->value += 10;
                } else if (comp_value != ref_value && ref_value <= 1) {
                    ++opp;
                }
            }
        }
        if (i > 1) {
            // Move the pointer up the tree (towards the root)
            comp_node = &(*(comp_node - comp_node->node_idx + comp_node->root_idx)); 
        }
    }
    bounds.same = same; bounds.opps = opp;
    return max((unsigned int)same, (unsigned int)opp) - min((unsigned int)same, (unsigned int)opp);
}

__device__ size_t* realloc(size_t& old_size, size_t new_size, size_t* old)
{
    size_t* new_array = (size_t*)malloc(new_size * sizeof(size_t));

    for (size_t i = 0; i < old_size; i++) new_array[i] = old[i];

    free(old);
    return new_array;
}

__device__ void add_alignments(const internal::Tree& tree, TreeNode* node, 
                               const size_t* aligned, const size_t last_unaligned_idx)
{
    size_t read_offset = 0, align_count = 0;
    size_t* alignments = new size_t[tree.snp_info[node->haplo_idx].end_index() - 
                                    tree.snp_info[node->haplo_idx].start_index()];
    size_t* values     = new size_t[tree.snp_info[node->haplo_idx].end_index() - 
                                    tree.snp_info[node->haplo_idx].start_index()];     

    for (size_t i = last_unaligned_idx; i < tree.reads; ++i) {
        if (i >= tree.snp_info[node->haplo_idx].start_index() &&
            i <= tree.snp_info[node->haplo_idx].end_index()  ) {
       
            // If the row crosses the snp
            if (tree.read_info[aligned[i]].start_index() <= node->haplo_idx &&
                tree.read_info[aligned[i]].end_index()   >= node->haplo_idx ) {
            
                read_offset  = tree.read_info[i].offset();
                auto element = tree.data[read_offset + node->haplo_idx];
                node->node_idx = aligned[i];
            
                if ((element == 0 && node->value == 0) || (element == 1 && node->value == 1)) {
                    values[align_count] = 1; alignments[align_count++] = aligned[i];
                    node->node_idx++;
                } else if ((element == 0 && node->value == 1) || (element == 1 && node->value == 0)) {
                    values[align_count] = 0; alignments[align_count++] = aligned[i];
                    node->node_idx++;
                }
            }
        }
    }

    // Move the found alignments to the node
    node->alignments     = align_count;
    node->read_ids       = new size_t[align_count];
    node->read_values    = new uint8_t[align_count];
    
    for (size_t i = 0; i < align_count; ++i) {
        node->read_ids[i] = alignments[i];
        node->read_values[i] = (uint8_t)values[i];
    }

    free(alignments); free(values);
}


__device__ void search_helper(const internal::Tree& tree, TreeNode* node, const size_t start_index, )
{
    // Set node parameters
    // Value only for DP
    node.haplo_idx = tree.search_snps[tree.last_searched_snp];
}

__global__ void search_tree(internal::Tree tree, size_t*  result, size_t* ubound)
{
    // Which alignments are and are not set
    size_t* set_alignments  = new size_t[tree.reads];
    size_t  last_aligned    = 0;
    unsigned int min_ubound = 0, unsigned int branches = 0;
    
    // initialize the alignments
    for (size_t i = 0; i < tree.reads; ++i) set_alignments[i] = i;

    // ------------------------------------- SEARCH START ---------------------------------------------------
    
    TreeNode node  = tree.node_manager.node(0);
    node.haplo_idx = tree.search_snps[tree.last_searched_snp];
    node.node_idx  = 0; node.value  = 0;
    node.lbound    = 0; node.ubound = *ubound;
    
    // Set the alignments for the tree root
    add_alignments(tree, &node, set_alignments, last_aligned);

    // Add the alignments to the overall alignments
    for (size_t i = last_aligned; i < last_aligned + node.alignments; ++i) {
            set_alignments[node.read_ids[i - last_aligned]] = set_alignments[i];
            set_alignments[i] = node.read_ids[i - last_aligned];
    } last_aligned += node.alignments;
    
    // Go over all the nodes that have not been searched and see how correlated they are
    // Parallel
    size_t temp = 0, res = 0;
    for (size_t i = tree.last_searched_snp + 1; i < tree.snps; ++i) {
        res = compare_snps(tree, &node, i);
        if (res > temp) *result = i;
    }

    // The first node has now been searched
    tree.last_searched_snp++;
    
    // Make the next 2 nodes point back to this one
    TreeNode left_child = tree.node_manager.node(1); TreeNode right_child = tree.node_manager.node(2);
    left_child.root_idx  = 0; right_child.root_idx  = 0;  
    left_child.value     = 0; right_child.value     = 1;
    left_child.node_idx  = 1; right_child.node_idx  = 2;
    left_child.haplo_idx = *result; right_child.haplo_idx = *result;
    
    // Now we can start the iterative search 

    // Do the alignments for the next nodes
    add_alignments(tree, &left_child , set_alignments, last_aligned);
    add_alignments(tree, &right_child, set_alignments, last_aligned);
    
    size_t a = compare_snps(tree, &left_child, 2);
    tree.last_searched_snp++;
    
//    *result = left_child.value;
 
    free(set_alignments);
}

}               // End namespace haplo
