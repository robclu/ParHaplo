#include "cuda.h"
#include "data.h"
#include "debug.h"
#include "graph_internal.h"
#include "math.h"

#define BLOCK_SIZE 1024

namespace haplo {

using graph_type = internal::Graph;
using data_type  = Data;

extern __shared__ size_t distance[];

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
    
    // Check if each of the values is valid
    if (read_info_one.start_index() <= snp_idx && read_info_one.end_index() >= snp_idx) {
        first_valid = true;
        first_value = data.data[read_info_one.offset() + snp_idx];
    }
    if (read_info_two.start_index() <= snp_idx && read_info_two.end_index() >= snp_idx) {
        second_valid = true;
        second_value = data.data[read_info_two.offset() + snp_idx];
    }
    
    // Load the distance into the shared array
    if (first_valid && second_valid) {
        if (first_value != second_value && first_value <= 1 && second_value <= 1) {
            distance[snp_idx] = 10;
            distance[snp_idx + threads] = 1;
        } else if (first_value != second_value && (first_value <= 1 || second_value <= 1)) {
            distance[snp_idx] = 5;
            distance[snp_idx + threads] = 1;
        } else if (first_value == second_value) {
            distance[snp_idx] = 0;
            distance[snp_idx + threads] = 1;
        }
    } else if (first_valid && !second_valid && first_value <= 1) {
        distance[snp_idx] = 5;
        distance[snp_idx + threads] = 1;
    } else if (second_valid && !first_valid && second_value <= 1) {
        distance[snp_idx] = 5;
        distance[snp_idx + threads] = 1;
    } else if (!first_valid && !second_valid) {
        distance[snp_idx] = 0;
        distance[snp_idx + threads] = 0;
    }
    __syncthreads();    
    
    // Set the values for the graph edge
    graph.edges[blockIdx.x].f1 = read_idx_one;
    graph.edges[blockIdx.x].f2 = read_idx_two;
    
    if (threadIdx.y == 0 && blockIdx.x == 2) {
        printf("R1: %i\n", read_idx_one);
        printf("R@: %i\n", read_idx_two);
        for (size_t i = 0; i < 6; ++i) 
            printf("%i ", distance[i]);
        printf("\n");
    } 
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
    graph.edges[blockIdx.x].distance = static_cast<float>(distance[0] / 10.f) /
                                       static_cast<float>(distance[threads])  ;
    
    if (threadIdx.y == 0 && blockIdx.x == 2) {
        for (size_t i = 0; i < 6; ++i) 
            printf("%i ", distance[i]);
        printf("\n");
    } 
}
 
__device__
void swap_edges(graph_type& graph, const size_t edge_idx_one, const size_t edge_idx_two)
{
    Edge temp = graph.edges[edge_idx_one];
    graph.edges[edge_idx_one] = graph.edges[edge_idx_two];
    graph.edges[edge_idx_two] = temp; 
}

// Sorts the edges
__device__ 
void bitonic_out_in_sort(graph_type& graph, const size_t block_idx, const size_t block_size)
{
    const size_t idx_one    = blockIdx.x + (block_idx * (block_size / 2));
    const size_t idx_two    = idx_one + (block_size - (blockIdx.x % (block_size / 2)) - 1) 
                            - (idx_one % (block_size / 2));
   
    // If the dimensions are in range
    if (idx_two < gridDim.x) {
        // The edges need to be swapped if the right one is larger than the left one
        // or if the left one has a value of 0.5, since those must be removed
        if ( (graph.edges[idx_one].distance < graph.edges[idx_two].distance                   ) ||
             (graph.edges[idx_one].distance == 0.5f && graph.edges[idx_two].distance != 0.5f) ) {
                swap_edges(graph, idx_one, idx_two);       
        }
    }
}

__device__
void bitonic_out_out_sort(graph_type& graph, const size_t block_idx, const size_t block_size)
{
    const size_t idx_one = blockIdx.y + (block_idx * (block_size / 2));
    const size_t idx_two = idx_one + (block_size / 2);

    // Check that the node index is in the first half of the block and the comp node is in range
    if (idx_two < gridDim.x) {
        if ( (graph.edges[idx_one].distance < graph.edges[idx_two].distance                   ) ||
             (graph.edges[idx_one].distance == 0.5f && graph.edges[idx_two].distance != 0.5f) ) {
                swap_edges(graph, idx_one, idx_two);       
        }        
    }    
}

// Reduces a level using a parallel bitonic sort
__device__
void reduce_edges(graph_type& graph)
{
    const size_t    passes     = static_cast<size_t>(ceil(log2(static_cast<double>(gridDim.x))));
    size_t          block_size = 2;
    
    if (blockIdx.x < gridDim.x / 2) {
        for (size_t pass = 0; pass < passes; ++pass) {
            // We need a logarithmically decreasing number of out-in passes 
            bitonic_out_in_sort(graph, blockIdx.x / (block_size / 2), block_size);
            //__syncthreads();

            // Now we need pass number of out-out bitonic sorts
            size_t out_out_block_size = block_size / 2; 
            for (size_t i = 0; i < pass; ++i) {
                bitonic_out_out_sort(graph, blockIdx.x / (out_out_block_size / 2), out_out_block_size);
                out_out_block_size /= 2;
                //__syncthreads();
            }
            block_size *= 2;
        }
    }
}

__device__ 
size_t find_last_valid_edge(const graph_type& graph, size_t last_edge_index)
{
    bool found_end = false;   
    while (!found_end) {
       if (graph.edges[last_edge_index].distance != 0.5f) 
            found_end = true;
        else 
            --last_edge_index;
    }
    return last_edge_index;
}

__global__
void search_graph(data_type data, graph_type graph, size_t threads) 
{
    const size_t total_vertices = data.reads * (data.reads - 1) / 2;
    
    map_distances(data, graph, threads);
    reduce_distances(data, graph, threads);
   
    reduce_edges(graph);
    
    const size_t last_valid_edge = find_last_valid_edge(graph, total_vertices - 1);
    
    if (threadIdx.y == 0 && blockIdx.x == 2) {
        for (size_t i = 0; i < data.reads * (data.reads -1) / 2; ++i)
            printf("%f ", graph.edges[i].distance);
        printf("\n");
        for (size_t i = 0; i < data.reads * (data.reads -1) / 2; ++i)
            printf("%i ", graph.edges[i].f1);
        printf("\n");        
        for (size_t i = 0; i < data.reads * (data.reads -1) / 2; ++i)
            printf("%i ", graph.edges[i].f2);
        printf("\n%i\n", last_valid_edge);    
    }
}

}               // End namespace haplo
