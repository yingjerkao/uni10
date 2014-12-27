UNI10_ROOT := /home/yundasie/uni10
UNI10_ROOT := /Users/yunda/myLibrary/uni10_dev
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a
HAMILTONIAN_INC := ../Hamiltonians

itebd.e: 2DiTEBD.cpp 2DiTEBD_tools.cpp
	$(CC) -I$(INC) -I$(HAMILTONIAN_INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<


clean:
	rm -f *.o *.e
