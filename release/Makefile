########################################################################################
#				 	                EXECUTABLE NAME 	                               #
########################################################################################

CXX_EXE 		:= parahaplo
NXX_EXE 		:= parahaplo

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
NXX_LDIR        := -L/usr/local/cuda-7.0/lib64

########################################################################################
#					                COMPILER FLAGS 					                   #
#                                                                                      #
# NOTE: To enable compiler warnings, remove -w and replace with :                      #
# 		--compiler-options -Wall                                                       #
########################################################################################

CXX_FLAGS 		:= -std=c++11 -Wall 
NXX_FLAGS       := -arch=sm_30 --std=c++11 -O3


########################################################################################
# 					                TARGET RULES 					                   #
#######################################################################################

.PHONY: all parahaplo evaluator build_parahaplo_and_run 

all: parahaplo
	
parahaplo: build_parahaplo
	
evaluator: build_evaluator
	
build_parahaplo_and_run: build_and_run

build_and_run: build_parahaplo
	./$(NXX_EXE)
	rm -rf *.o

evaluator.o: ../haplo/evaluator.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<
	
evaluator_main.o: evaluator_main.cpp 
	$(CXX) $(CXX_INCLUDE) $(CXX_FLAGS) -o $@ -c $<

parahaplo.o: parahaplo.cu
	$(NXX) $(NXX_INCLUDE) $(NXX_FLAGS) -o $@ -dc $<

build_evaluator: evaluator.o evaluator_main.o 
	$(CXX) -o $(CXX_EXE) $+ $(CXX_LDIR) $(CXX_LIBS)	
	
build_parahaplo: NXX_FLAGS += -DSTAND_ALONE
build_parahaplo: parahaplo.o 
	$(NXX) -o $(NXX_EXE) $+ $(NXX_LDIR) $(NXX_LIBS)	

clean:
	rm -rf *.o
	rm -rf $(CXX_EXE) 
	rm -rf $(NXX_EXE) 

