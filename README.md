*A parallel implementation of the haplotype assembly problem*

# Parahaplo

Provides a parallel implementation of the haplotype assembly problem, that works with a multiple parallel
hardware:

* __GPU__: This is the intended usage and provides the best performance, it uses CUDA.
* __CPU__: Parallel CPU implementation using Intel TBB.
* __PHI__: Intel Xeon Phi implementation, it uses Intel TBB

#Organisation 

The repository is structured as follows:

```
|--- parahaplo
|    |--- LICENCE.md       
|    |--- README.md
|    |--- doc/                                  (documentation)
|    |    |--- prelim/                          (preliminary report)
|    |    |--- final/                           (final report)
|    |--- haplo/                                (source code)
|    |--- ref/                                  (reference material)
|    |    |--- links.md                         (links to any relevant information)
|    |--- tests/                                (unit tests for library)
|    |    |--- block_tests.cpp                  (block tests)
|    |    |--- equality_checker_tests.cpp       (tests for checking if rows/cols are equal)
|    |    |--- Makefile                         (Makefile for tests)
|    |    |--- small_container_tests.cpp        (tests for small (binary) containers)
|    |    |--- tests.cpp                        (main test file for running all tests)
|    |    |--- unsplittable_block_tests.cpp     (tests for unsplittable blocks)
|    |    |--- input_files/                     (input files for tests)
```

# Dependencies 

# Compiling 
