SRC := ../src2/
INC := ../src2/
LIB := ../src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

all:
	rm -f egU1.e
	make egU1.e

egU1.e: egU1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egB1.e: egB1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
scl.e: Scalar.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
trc.e: Trace.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)




clean:
	rm -f *.o *.e asd
