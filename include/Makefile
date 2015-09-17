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

CXX_FLAGS 		:= -std=c++11 -w -g
NXX_FLAGS       := -arch=sm_30
PXX_FLAGS       :=

########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: tests clean asm

asm: CXX_FLAGS += -S -fverbose-asm 
asm: make_asm
	
parahaplo_tests.s: parahaplo_tests.cpp
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
make_asm: parahaplo_tests.s
	as -alhnd parahaplo_tests.s > parahaplo_tests.lst 
	
tests: build_and_run

build_and_run: build_tests
	rm -rf *.o
	./$(CXX_EXE) --log_level=test_suite
	
build: build_tests
	
parahaplo_tests.o: parahaplo_tests.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

build_tests: parahaplo_tests.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)
	
clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 

clean_asm:
	rm -rf parahaplo_tests.s
	rm -rf parahaplo_tests.lst
	rm -rf a.out

