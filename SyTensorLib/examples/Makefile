SRC := ../src2/
INC := ../src2/
LIB := ../src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

all:
	rm -f egQ1.e
	make egQ1.e

egQ1.e: egQ1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
scl.e: Scalar.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
trc.e: Trace.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)




clean:
	rm -f *.o *.e asd
