########################################################################################
#				 	                EXECUTABLE NAME 	                               #
########################################################################################

CXX_EXE 		:= parahaplo_tests

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

CXX 			:= clang++ 
NXX             := nvcc --ccbin $(CXX)

########################################################################################
#				                    INCLUDE DIRECTORIES 		                       #
########################################################################################

CXX_INCLUDE     := -I$(TBBROOT)/include -I/Users/sashatn/Documents/boost_1_59_0/include

########################################################################################
#				   	                   LIBRARIES 						               #
########################################################################################

CXX_LIBS 		:= -lboost_unit_test_framework -lboost_iostreams -lm -ltbb
CXX_LDIR  	    := 

########################################################################################
#					                COMPILER FLAGS 					                   #
#                                                                                      #
# NOTE: To enable compiler warnings, remove -w and replace with :                      #
# 		--compiler-options -Wall                                                       #
########################################################################################

CXX_FLAGS 		:= -std=c++11 -Wall -g
NXX_FLAGS       := -arch=sm_30
PXX_FLAGS       :=

ALL_TESTS       =   small_container_tests.o             \
					equality_checker_tests.o            \
					data_converter.o                    \
					data_converter_tests.o              \
					block_tests.o                       \
					tree_tests.o                        \
					node_container_tests.o              \
					unsplittable_block_tests.o          \
					tests.o 

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: all block_tests container_tests checker_tests converter_tests node_tests tree_tests \
		unsplittable_block_tests clean

all: build_and_run

build_and_run: build_tests
	./$(CXX_EXE) --log_level=test_suite
	rm -rf *.o
	
build: build_tests

block_tests.o: block_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

equality_checker_tests.o: equality_checker_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

data_converter.o: ../haplo/data_converter.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
data_converter_tests.o: data_converter_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

node_container_tests.o: node_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
small_container_tests.o: small_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

tree_tests.o: tree_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

unsplittable_block_tests.o: unsplittable_block_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

tests.o: tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

build_tests: $(ALL_TESTS)
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)

container_tests: CXX_FLAGS += -DSTAND_ALONE
container_tests: small_container_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)
	
block_tests: CXX_FLAGS += -DSTAND_ALONE
block_tests: block_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	
	
checker_tests: CXX_FLAGS += -DSTAND_ALONE
checker_tests: equality_checker_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

converter_tests: CXX_FLAGS += -DSTAND_ALONE
converter_tests: data_converter.o data_converter_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

node_tests: CXX_FLAGS += -DSTAND_ALONE
node_tests: node_container_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

tree_tests: CXX_FLAGS += -DSTAND_ALONE
tree_tests: tree_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	
	
unsplittable_block_tests: CXX_FLAGS += -DSTAND_ALONE
unsplittable_block_tests: unsplittable_block_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	

clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 

