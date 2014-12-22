UNI10_ROOT := /Users/yunda/myLibrary/uni10_dev
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
HAMILTONIAN_INC := ../Hamiltonians

itebd.e: iTEBD.cpp iTEBD_tools.cpp
	$(CC) -I$(INC) -I$(HAMILTONIAN_INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<

clean:
	rm -f *.o *.e
