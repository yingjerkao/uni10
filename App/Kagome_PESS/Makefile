#UNI10_ROOT := /home/yundasie/uni10
UNI10_ROOT := /Users/Mac/guni10Lib
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a
KagomeLib := KagomePESSLib

KagomePESS: KagomePESS.cpp
	$(CC) -I$(KagomeLib) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
Measuretest: Measuretest.cpp
	$(CC) -I$(KagomeLib) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
Hamiltonian: Hamiltonian.cpp
	$(CC) -I$(KagomeLib) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<

clean:
	rm -f *.o *.e KagomePESS Hamiltonian Measuretest
