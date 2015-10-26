#include "cuda.h"
#include "data.h"
#include "debug.h"
#include "graph_internal.h"
#include "math.h"

#define BLOCK_SIZE 1024
//#undef DEBUG

namespace haplo {

using graph_type = internal::Graph;
using data_type  = Data;

extern __shared__ size_t distance[];
extern __shared__ uint8_t shared_sets[];

__device__
void map_distances(data_type& data, graph_type& graph, const size_t threads)
{
    const size_t total_elems  = data.reads * (data.reads - 1) / 2;
    const size_t read_idx_one = blockIdx.x / (data.reads - 1);
    const size_t subt_elems   = (data.reads - read_idx_one) * (data.reads - read_idx_one - 1) / 2;
    const size_t read_idx_two = read_idx_one + (blockIdx.x - (total_elems - subt_elems)) + 1;
    const size_t snp_idx      = threadIdx.y;

    const auto read_info_one = data.read_info[read_idx_one];    
    const auto read_info_two = data.read_info[read_idx_two];    
    
    bool    first_valid = false , second_valid = false;
    uint8_t first_value = 0     , second_value = 0    ;

    if (read_idx_one < data.reads && read_idx_two < data.reads) {    
        // Check if each of the values is valid
        if (read_info_one.start_index() <= snp_idx && read_info_one.end_index() >= snp_idx) {
            first_valid = true;
            first_value = data.data[read_info_one.offset() + snp_idx - read_info_one.start_index()];
        }
        if (read_info_two.start_index() <= snp_idx && read_info_two.end_index() >= snp_idx) {
            second_valid = true;
            second_value = data.data[read_info_two.offset() + snp_idx - read_info_two.start_index()];
        }
        
        // Set the values for the graph edge
        graph.edges[blockIdx.x].f1 = read_idx_one;
        graph.edges[blockIdx.x].f2 = read_idx_two;
    } 
    
    // Load the distance into the shared array
    if (first_valid && second_valid) {
        if (first_value != second_value && first_value <= 1 && second_value <= 1) {
            distance[snp_idx] = 10; distance[snp_idx + threads] = 1;
        } else if (first_value != second_value && (first_value <= 1 || second_value <= 1)) {
            distance[snp_idx] = 5; distance[snp_idx + threads] = 1;
        } else if (first_value == second_value && first_value <= 1 && second_value <= 1) {
            distance[snp_idx] = 0; distance[snp_idx + threads] = 1;
        } else if (first_value == second_value & first_value > 1 && second_value > 1) {
            distance[snp_idx] = 0; distance[snp_idx + threads] = 0;
        }   
    } else if (first_valid && !second_valid && first_value <= 1) {
        distance[snp_idx] = 5; distance[snp_idx + threads] = 1;
    } else if (second_valid && !first_valid && second_value <= 1) {
        distance[snp_idx] = 5; distance[snp_idx + threads] = 1;
    } else if (!first_valid && !second_valid) {
        distance[snp_idx] = 0; distance[snp_idx + threads] = 0;
    }
    __syncthreads();    
}

__device__
void reduce_distances(data_type& data, graph_type& graph, const size_t threads) 
{
    size_t reduction_threads = threads;
    
    while (reduction_threads > 1) {
        size_t other_idx_one = threadIdx.y + (reduction_threads / 2);
        size_t other_idx_two = other_idx_one + threads;
        
        if (threadIdx.y < (reduction_threads / 2)) {                
            distance[threadIdx.y] += distance[other_idx_one];
            distance[threadIdx.y + threads] += distance[other_idx_two];
        }
        if (reduction_threads % 2 == 1) {
            // Load the extra value
            if (threadIdx.y == reduction_threads / 2) {
                distance[threadIdx.y] = distance[threadIdx.y + reduction_threads / 2];
                distance[threadIdx.y + threads] = distance[threadIdx.y + threads + reduction_threads / 2];
            }
            reduction_threads /= 2; reduction_threads += 1;
        } else reduction_threads /= 2;
        __syncthreads();
    }
    // Set the weight of the edge
    if (distance[threads] > 0 && threadIdx.y == 0) {
        graph.edges[blockIdx.x].distance = static_cast<float>(distance[0] / 10.f) /
                                           static_cast<float>(distance[threads])  + 0.5f;
    } else if (threadIdx.y == 0 && distance[threads] == 0) {
        graph.edges[blockIdx.x].distance = 1.0f;
    }
    
    // Modify the value so that it's easier to sort 
    if (graph.edges[blockIdx.x].distance == 1.0f && threadIdx.y == 0) {
        graph.edges[blockIdx.x].distance = 0.0f;
    }
}
 
__device__
void swap_edges(graph_type& graph, const size_t edge_idx_one, const size_t edge_idx_two)
{
    Edge temp = graph.edges[edge_idx_one];
    graph.edges[edge_idx_one] = graph.edges[edge_idx_two];
    graph.edges[edge_idx_two] = temp; 
}

// Out to in operation for bitonic sort
__global__ 
void bitonic_out_in_sort(graph_type graph, const size_t block_size, const size_t total_elements)
{
    const size_t block_idx  = blockIdx.x / (block_size / 2);
    const size_t idx_one    = blockIdx.x + (block_idx * (block_size / 2));
    const size_t idx_two    = idx_one + (block_size - (blockIdx.x % (block_size / 2)) - 1) 
                            - (idx_one % (block_size / 2));
    
    // If the dimensions are in range
    if (idx_one < total_elements && idx_two < total_elements && threadIdx.y == 0) {
        // The edges need to be swapped if the right one is larger than the left one
        // or if the left one has a value of 0.5, since those must be removed
        if (graph.edges[idx_one].distance <= graph.edges[idx_two].distance) {
            swap_edges(graph, idx_one, idx_two);
        } 
    }
}

// Out to out operation fr bitonic sort
__global__
void bitonic_out_out_sort(graph_type graph, const size_t block_size, const size_t total_elements)
{
    const size_t block_idx = blockIdx.x / (block_size / 2);
    const size_t idx_one   = blockIdx.x + (block_idx * (block_size / 2));
    const size_t idx_two   = idx_one + (block_size / 2);

    // Check that the node index is in the first half of the block and the comp node is in range
    if (idx_one < total_elements && idx_two < total_elements && threadIdx.y == 0) {
        // The edges need to be swapped if the right one is larger than the left one
        // or if the left one has a value of 0.5, since those must be removed
        if (graph.edges[idx_one].distance <= graph.edges[idx_two].distance) {
            swap_edges(graph, idx_one, idx_two);
        }
    }    
}

__device__ 
size_t find_last_valid_edge(const graph_type& graph)
{
    bool   found_end       = false;   
    size_t last_valid_edge = 0;
    while (!found_end) {
       if (graph.edges[last_valid_edge].distance != 0.0f) 
            ++last_valid_edge;
        else 
            found_end = true;
    }
    return last_valid_edge;
}

// Prints the edges of a node
__global__
void print_edges(data_type data, graph_type graph)
{
    if (threadIdx.y == 0 && blockIdx.x == 0) {   
        for (size_t i = 0; i < data.reads * (data.reads -1) / 2; ++i) {
            if (graph.edges[i].distance != 0.0f) {
                printf("%.4f ", graph.edges[i].distance);
                printf("%i ", graph.edges[i].f1);
                printf("%i\n", graph.edges[i].f2);
            }
        }
    } 
}

__global__
void search_graph(data_type data, graph_type graph, size_t threads) 
{
    const size_t total_vertices = data.reads * (data.reads - 1) / 2;
    
    map_distances(data, graph, threads);
    reduce_distances(data, graph, threads);
}

// For now we assume the number of partitions is less than 6000
// Add dynamic parallelism
template <uint8_t Set> __device__ 
uint8_t in_set(graph_type&  graph        , const uint32_t fragment, const size_t block_idx_start, 
               const size_t block_idx_end, uint8_t*      shared_set)
{ 
    // Need to extend this to use multiple blocks 
    const size_t frag_idx           = threadIdx.y;
    size_t       reduction_threads  = Set == 1 ? graph.set_one_size : graph.set_two_size;

    if (blockIdx.x == 0 && threadIdx.y == 0) {
        for (size_t i = 0; i < reduction_threads; ++i) {
            if (Set == 1 ? graph.set_one[i] == fragment : graph.set_two[i] == fragment) {
                return true;
            }
        }
        return false;
    }
/*
    if (blockIdx.x >= block_idx_start && blockIdx.x < block_idx_end && frag_idx < reduction_threads) {    
        shared_set[frag_idx] = Set == 1 ? graph.set_one[frag_idx] == fragment 
                                        : graph.set_two[frag_idx] == fragment;
        __syncthreads();
        
        // Perform the reduction  
        while (reduction_threads > 1) {
            if (frag_idx < (reduction_threads / 2)) {
                // If the thread is part of the first half of the threads
                shared_set[frag_idx] += shared_set[frag_idx + (reduction_threads / 2)];
            } 
            if (reduction_threads % 2 == 1) {
                if (frag_idx == reduction_threads / 2) {
                    // If there are an odd number of elements then move the last element
                    shared_set[frag_idx] = shared_set[frag_idx + (reduction_threads / 2)];
                }
                reduction_threads /= 2; reduction_threads += 1;
            } else reduction_threads /= 2;               
            __syncthreads();
        }
    }
    if (threadIdx.y == 0 && blockIdx.x == 0) result = shared_set[0];
    if (frag_idx < Set == 1 ? graph.set_one_size : graph.set_two_size) shared_set[frag_idx] = 0;
    __syncthreads();
    
    return result;
*/
}



__device__
void partition_next_largest_fragment(graph_type& graph, size_t& last_set_edge, const size_t last_valid_edge,
                                     uint8_t* shared_set)
{
    const size_t initial_edge = last_set_edge;
    bool   found              = false;
   
    if (blockIdx.x == 0 && threadIdx.y == 0) { 
        while (last_set_edge < last_valid_edge - 1 && !found) {
            uint8_t f1_in_set_1 = in_set<1>(graph, graph.edges[last_set_edge].f1, 0, 1, shared_set); 
            uint8_t f1_in_set_2 = in_set<2>(graph, graph.edges[last_set_edge].f1, 0, 1, shared_set); 
            uint8_t f2_in_set_1 = in_set<1>(graph, graph.edges[last_set_edge].f2, 0, 1, shared_set); 
            uint8_t f2_in_set_2 = in_set<2>(graph, graph.edges[last_set_edge].f2, 0, 1, shared_set);  
            
            if (f1_in_set_1 && !f2_in_set_2 && !f2_in_set_1) {
                // f1 in set one, add f2 to set two
                graph.set_two[graph.set_two_size] = graph.edges[last_set_edge].f2;
                ++graph.set_two_size;
                found = true;
            } else if (f2_in_set_1 && !f1_in_set_2 && !f1_in_set_1) {
                // f2 in set one, add f1 to set two
                graph.set_two[graph.set_two_size] = graph.edges[last_set_edge].f1;
                ++graph.set_two_size;
                found = true;
            } else if (f1_in_set_2 && !f2_in_set_1 && !f2_in_set_2) {
                // f1 in set 2, add f2 so set one
                graph.set_one[graph.set_one_size] = graph.edges[last_set_edge].f2;
                ++graph.set_one_size;
                found = true;
            } else if (f2_in_set_2 && !f1_in_set_1 && !f1_in_set_2) {
                // f2 in set 2, add f1 to set one
                graph.set_one[graph.set_one_size] = graph.edges[last_set_edge].f1;
                ++graph.set_one_size;
                found = true;
            } else ++last_set_edge;
        }
        if (found) {
            Edge temp                   = graph.edges[initial_edge];
            graph.edges[initial_edge]   = graph.edges[last_set_edge];
            graph.edges[last_set_edge]  = temp;
            last_set_edge               = initial_edge + 1;
        } else last_set_edge = initial_edge;
    }
}

__device__
void partition_next_smallest_fragment(graph_type& graph, size_t& last_set_edge, uint8_t* shared_set)
{
    const size_t initial_edge = last_set_edge;
    bool   found = false;
    
    if (blockIdx.x == 0 && threadIdx.y == 0) { 
        while (last_set_edge >= 1 && !found) {
            uint8_t f1_in_set_1 = in_set<1>(graph, graph.edges[last_set_edge].f1, 0, 1, shared_set); 
            uint8_t f1_in_set_2 = in_set<2>(graph, graph.edges[last_set_edge].f1, 0, 1, shared_set); 
            uint8_t f2_in_set_1 = in_set<1>(graph, graph.edges[last_set_edge].f2, 0, 1, shared_set); 
            uint8_t f2_in_set_2 = in_set<2>(graph, graph.edges[last_set_edge].f2, 0, 1, shared_set); 
            
            if (f1_in_set_1 && !f2_in_set_1 && !f2_in_set_2) {
                // f1 in set one, add f2 to set one
                graph.set_one[graph.set_one_size] = graph.edges[last_set_edge].f2;
                ++graph.set_one_size;
                found = true;
            } else if (f2_in_set_1 && !f1_in_set_1 && !f1_in_set_2) {
                // f2 in set one, add f1 to set one
                graph.set_one[graph.set_one_size] = graph.edges[last_set_edge].f1;
                ++graph.set_one_size;
                found = true;
            } else if (f1_in_set_2 && !f2_in_set_2 && !f2_in_set_1) {
                // f1 in set 2, add f2 so set two
                graph.set_two[graph.set_two_size] = graph.edges[last_set_edge].f2;
                ++graph.set_two_size;
                found = true;
            } else if (f2_in_set_2 && !f1_in_set_2 && !f1_in_set_1) {
                // f2 in set 2, add f1 to set one
                graph.set_two[graph.set_two_size] = graph.edges[last_set_edge].f1;
                ++graph.set_two_size;
                found = true;
            } else --last_set_edge;
        }
        if (found) {
            Edge temp                   = graph.edges[initial_edge];
            graph.edges[initial_edge]   = graph.edges[last_set_edge];
            graph.edges[last_set_edge]  = temp;
            last_set_edge               = initial_edge - 1; 
        } else last_set_edge = initial_edge;
    }
}

__device__ 
void print_sets(graph_type& graph)
{
    if (threadIdx.y == 0 && blockIdx.x == 0) {
        for (size_t i = 0; i < graph.set_one_size; ++i) 
            printf("%i ", graph.set_one[i]);
        printf("\n");
        for (size_t i = 0; i < graph.set_two_size; ++i) 
            printf("%i ", graph.set_two[i]);
        printf("\n"); 
    }
}

__global__ 
void map_to_partitions(data_type data, graph_type graph)
{
    const size_t last_valid_edge        = find_last_valid_edge(graph);
    size_t       last_set_edge_forward  = 1;
    size_t       last_set_edge_backward = last_valid_edge - 1;
    bool         keep_partitioning      = true;
   
    extern __shared__ uint8_t sets[];
    
#ifdef DEBUG
    if (threadIdx.y == 0 && blockIdx.x == 0)
        printf("LVE : %i\n", last_valid_edge);
#endif
    
    // Add the first elements in the partitions
    graph.set_one[0] = graph.edges[0].f1 + 1; graph.set_one_size = 1;
    graph.set_two[0] = graph.edges[0].f2 + 1; graph.set_two_size = 1;
    
    // Partition the remaining fragments
    if (threadIdx.y == 0 && blockIdx.x == 0) {
        while (keep_partitioning) {
            size_t last_edge_back_before = last_set_edge_backward;
            size_t last_edge_fwd_before  = last_set_edge_forward;

            partition_next_largest_fragment(graph, last_set_edge_forward, last_valid_edge, &sets[0]);
            partition_next_smallest_fragment(graph, last_set_edge_backward, &sets[0]);   
            
            if (last_set_edge_forward == last_valid_edge || last_set_edge_backward == 1                           ||
                (last_set_edge_forward == last_edge_fwd_before && last_set_edge_backward == last_edge_back_before) ||
                (graph.set_one_size + graph.set_two_size == data.reads)                                           ) {
                keep_partitioning = false;
            }
        }  
    }
    
    print_sets(graph);
}

}               // End namespace haplo
