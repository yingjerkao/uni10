UNI10_ROOT := /Users/yunda/myLibrary/uni10_dev
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
OBJ := ham.o

main.e: main.cpp ham.o
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $< $(OBJ)
ham.o: hamiltonian.cpp
	$(CC) -c -I$(INC) $(FLAGS) $< -o $@

clean:
	rm -f *.o *.e
