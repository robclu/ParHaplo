########################################################################################
#				 	                EXECUTABLE NAME 	                               #
########################################################################################

CXX_EXE 		:= parahaplo_test_suite

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
NXX             := nvcc --ccbin $(CXX)

########################################################################################
#				                    INCLUDE DIRECTORIES 		                       #
########################################################################################

CXX_INCLUDE     := -I$(TBBROOT)/include

########################################################################################
#				   	                   LIBRARIES 						               #
########################################################################################

CXX_LIBS 		:= -lboost_unit_test_framework -lboost_iostreams -lm \
					-Wl,-rpath,`$(TBB_LIBRARY_RELEASE) -L$(TBB_LIBRARY_RELEASE)` -ltbb	
CXX_LDIR  	    := 

########################################################################################
#					                COMPILER FLAGS 					                   #
#                                                                                      #
# NOTE: To enable compiler warnings, remove -w and replace with :                      #
# 		--compiler-options -Wall                                                       #
########################################################################################

CXX_FLAGS 		:= -std=c++11 -Wall -O3
NXX_FLAGS       := -arch=sm_30
PXX_FLAGS       :=

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: all tests clean

all: tests
	
tests: build_and_run

build_and_run: build_tests
	rm -rf *.o
	./$(CXX_EXE) --log_level=test_suite
	
build: build_tests

small_container_tests.o: small_container_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

block_tests.o: block_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

build_tests: small_container_tests.o block_tests.o tests.o
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)
	
clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 

