# paraHaplo
paraHaplo is an implementation of the parallel individual haplotyping assembly problem. 

#Description 

This is for my final year project, which is an investigation into the parallelization of the individual haplotyping assembly problem. There will be 3 different implementations:
* __Multi-threaded C++__: Using standard C++ and the SIMD instructions.
* __CUDA__: Parallel implementation using CUDA, capable of running on Nvidia graphics cards. 
* __Xeon Phi__: Parallel implementation capable of running on the Xeon Phi.

## Languages

C++ (11) is used for all implementations. The Nvidia CUDA API is used for the CUDA implementation

#Organisation 

The repository is structured as follows:

```
| ROOT
-| doc (documentation)
--| prelim (preliminary report)
--| final (final report)
-| src (source code)
-| ref (reference material)
--| links (links to any relevant information)
```
