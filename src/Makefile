########################################################################################
#				 	                EXECUTABLE NAME 	                               #
########################################################################################

CXX_EXE 		:= parahaplo_cpp_test_suite

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
#	PXX     : The compiler for the intel xeon phi implementation                       #
########################################################################################

CXX 			:= g++ 
NXX             := nvcc --ccbin $(CXX)
PXX             :=

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

CXX_FLAGS 		:= -std=c++11 -w -g
NXX_FLAGS       := -arch=sm_30
PXX_FLAGS       :=

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: cpp clean_cpp asm

asm: CXX_FLAGS += -S -fverbose-asm 
asm: make_asm
	
parahaplo_cpp_tests.s: parahaplo_cpp_tests.cpp
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
make_asm: parahaplo_cpp_tests.s
	as -alhnd parahaplo_cpp_tests.s > parahaplo_cpp_tests.lst 
	
cpp: build_and_run_cpp

build_and_run_cpp: build_cpp_tests
	rm -rf *.o
	./$(CXX_EXE) --log_level=test_suite
	
build: build_tests
	
parahaplo_cpp_tests.o: parahaplo_cpp_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

build_cpp_tests: parahaplo_cpp_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)
	
clean_cpp:
	rm -rf *.o
	rm -rf $(CXX_EXE) 

clean_asm:
	rm -rf parahaplo_cpp_tests.s
	rm -rf parahaplo_cpp_tests.lst
	rm -rf a.out

