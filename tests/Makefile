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

CXX 			:= g++-4.8
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
NXX_FLAGS       := -arch=sm_30 --std=c++11 -O3
PXX_FLAGS       :=

ALL_TESTS       =   small_container_tests.o             \
					data_converter.o                    \
					data_converter_tests.o              \
					evaluator.o                         \
					evaluator_tests.o                   \
					block_tests.o                       \
					subblock_tests.o                    \
					tests.o 

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: all

all: build_and_run

build_and_run: build_tests
	./$(CXX_EXE) --log_level=test_suite
	rm -rf *.o
	
build: build_tests

block_tests.o: block_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

data_converter.o: ../haplo/data_converter.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
data_converter_tests.o: data_converter_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

evaluator.o: ../haplo/evaluator.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
evaluator_tests.o: evaluator_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
small_container_tests.o: small_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

subblock_tests.o: subblock_tests.cpp
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

tests.o: tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

block_tests: CXX_FLAGS += -DSTAND_ALONE
block_tests: block_tests.o 
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

subblock_tests: CXX_FLAGS += -DSTAND_ALONE
subblock_tests: subblock_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	
	
snp_info_tests.o: snp_info_tests.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -c $<

snp_info_tests: NXX_FLAGS += -DSTAND_ALONE
snp_info_tests: snp_info_tests.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

build_tests: $(ALL_TESTS)
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)

clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 
	rm -rf $(NXX_EXE) 

