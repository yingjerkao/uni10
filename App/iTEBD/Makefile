UNI10_ROOT := /home/yundasie/uni10
#UNI10_ROOT := /Users/yunda/Research/TensorNetwork/uni10
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a

itebd_u1f.e: iTEBD_U1_fixed.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
titebd_u1f.e: TiTEBD_U1_fixed.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
itebd.e: iTEBD.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<

clean:
	rm -f *.o *.e asd
