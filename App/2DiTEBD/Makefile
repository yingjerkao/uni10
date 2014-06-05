UNI10_ROOT := /home/yundasie/uni10
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
OBJ := $(LIB)/libuni10.a

itebd.e: 2DiTEBD.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd
