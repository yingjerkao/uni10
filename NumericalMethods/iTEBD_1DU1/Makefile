UNI10_ROOT := /Users/yunda/myLibrary/uni10_dev
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC  := g++
FLAGS:= -O3 -std=c++11
LIBS += -lblas -llapack -lm -luni10
HAMILTONIAN_INC := ../Hamiltonians

itebd.e: iTEBD_U1.cpp iTEBD_U1_tools.cpp
	$(CC) $(FLAGS) -I$(INC) -I$(HAMILTONIAN_INC) -L$(LIB) $(LIBS) $< -o $@


clean:
	rm -f *.o *.e
