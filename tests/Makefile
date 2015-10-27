########################################################################################
#				 	                EXECUTABLE NAME 	                               #
########################################################################################

CXX_EXE 		:= parahaplo_tests
NXX_EXE 		:= parahaplo_cuda

########################################################################################
#                                   INSTALL DIRECTORY                                  #
#                                                                                      #
# Note: Will implement this later, if necessary - which I doubt                        #
########################################################################################

INSTALL_DIR     := /usr/include
UNINSTALL_DIR   := $(addsuffix /nano, $(INSTALL_DIR))

########################################################################################
#					                  COMPILERS						                   #
#					                                                                   #
#	CXX     : The compiler for the c++ implementation                                  #
#	NXX     : The compiler for the nvidia cuda gpu implementation                      #
########################################################################################

CXX 			:= g++
NXX             := nvcc -ccbin $(CXX)

########################################################################################
#				                    INCLUDE DIRECTORIES 		                       #
########################################################################################

CXX_INCLUDE     := -I$(TBBROOT)/include -I/Users/sashatn/Documents/boost_1_59_0/include
NXX_INCLUDE     := -I/usr/local/cuda-7.0/include
CXX_INCLUDE     += $(NXX_INCLUDE)

########################################################################################
#				   	                   LIBRARIES 						               #
########################################################################################

CXX_LIBS 		:= -lboost_unit_test_framework -lboost_iostreams -lm -ltbb
CXX_LDIR  	    := -L/Users/sashatn/Documents/boost_1_59_0/lib
NXX_LIBS        := -lcudart -lcuda -lcudadevrt $(CXX_LIBS)
NXX_LDIR        := -L/usr\local/cuda-7.0/lib64

########################################################################################
#					                COMPILER FLAGS 					                   #
#                                                                                      #
# NOTE: To enable compiler warnings, remove -w and replace with :                      #
# 		--compiler-options -Wall                                                       #
########################################################################################

CXX_FLAGS 		:= -std=c++11 -Wall 
NXX_FLAGS       := -arch=sm_30 --std=c++11 -G -g
PXX_FLAGS       :=

ALL_TESTS       =   small_container_tests.o             \
					comparator_tests.o                  \
					equality_checker_tests.o            \
					data_converter.o                    \
					data_converter_tests.o              \
					evaluator_test.o                    \
					graph_tests.o                       \
					block_tests.o                       \
					sorter_tests.o                      \
					node_selector_tests.o               \
					search_node_tests.o                 \
					node_manager_tests.o                \
					tree_tests.o                        \
					node_container_tests.o              \
					link_container_tests.o              \
					bounder_tests.o                     \
					subblock_tests.o                    \
					tests.o 

CUDA_TESTS      = cuda_container_tests.o                

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: all block_tests container_tests checker_tests converter_tests node_constiner_tests      \
		node_selector_tests link_tests tree_tests comparator_tests sorter_tests                 \
		search_node_tests node_manager_tests subblock_tests evaluator_tests                     \
		cuda_container_tests clean

all: build_and_run

build_and_run: build_tests
	./$(CXX_EXE) --log_level=test_suite
	rm -rf *.o
	
build: build_tests

block_tests.o: block_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

bounder_tests.o: bounder_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
comparator_tests.o: comparator_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

equality_checker_tests.o: equality_checker_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

data_converter.o: ../haplo/data_converter.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
data_converter_tests.o: data_converter_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

evaluator.o: ../haplo/evaluator.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
evaluator_tests.o: evaluator_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
link_container_tests.o: link_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

node_container_tests.o: node_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

node_selector_tests.o: node_selector_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

search_node_tests.o: search_node_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

node_manager_tests.o: node_manager_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
small_container_tests.o: small_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

sorter_tests.o: sorter_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

tree_tests.o: tree_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

subblock_tests.o: subblock_tests.cpp
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

tests.o: tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

block_tests: CXX_FLAGS += -DSTAND_ALONE
block_tests: block_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

bounder_tests: CXX_FLAGS += -DSTAND_ALONE
bounder_tests: bounder_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

checker_tests: CXX_FLAGS += -DSTAND_ALONE
checker_tests: equality_checker_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

comparator_tests: CXX_FLAGS += -DSTAND_ALONE
comparator_tests: comparator_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)

container_tests: CXX_FLAGS += -DSTAND_ALONE
container_tests: small_container_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)

converter_tests: CXX_FLAGS += -DSTAND_ALONE
converter_tests: data_converter.o data_converter_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	
	
evaluator_tests: CXX_FLAGS += -DSTAND_ALONE
evaluator_tests: evaluator.o evaluator_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

link_tests: CXX_FLAGS += -DSTAND_ALONE
link_tests: link_container_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

node_container_tests: CXX_FLAGS += -DSTAND_ALONE
node_container_tests: node_container_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

node_selector_tests: CXX_FLAGS += -DSTAND_ALONE
node_selector_tests: node_selector_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

search_node_tests: CXX_FLAGS += -DSTAND_ALONE
search_node_tests: search_node_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

node_manager_tests: CXX_FLAGS += -DSTAND_ALONE
node_manager_tests: node_manager_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

sorter_tests: CXX_FLAGS += -DSTAND_ALONE
sorter_tests: sorter_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

subblock_tests: CXX_FLAGS += -DSTAND_ALONE
subblock_tests: subblock_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

tree_tests: CXX_FLAGS += -DSTAND_ALONE
tree_tests: tree_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

build_tests: $(ALL_TESTS)
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)

cuda_container_tests.o: cuda_container_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -c $<

cuda_graph_tests.o: cuda_graph_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -dc $<
	
cuda_link_tests.o: cuda_link_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -c $<

cuda_tree_tests.o: cuda_tree_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -dc $<

snp_info_tests.o: snp_info_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -c $<

tree_kernels_gpu.o: ../haplo/tree_kernels_gpu.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -dc $<

cuda_container_tests: NXX_FLAGS += -DSTAND_ALONE
cuda_container_tests: cuda_container_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

cuda_graph_tests: NXX_FLAGS += -DSTAND_ALONE
cuda_graph_tests: cuda_graph_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	
	
cuda_link_tests: NXX_FLAGS += -DSTAND_ALONE
cuda_link_tests: cuda_link_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

cuda_tree_tests: NXX_FLAGS += -DSTAND_ALONE
cuda_tree_tests: cuda_tree_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

snp_info_tests: NXX_FLAGS += -DSTAND_ALONE
snp_info_tests: snp_info_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 
	rm -rf $(NXX_EXE) 

