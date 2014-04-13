SRC := ../src2/
INC := /opt/uni10/Library/uni10-0.1.0/include/
LIB := /opt/uni10/Library/uni10-0.1.0/lib/
#OBJ := $(LIB)libtensorF.a
#LINK := -mkl=parallel -luni10
OPT:=-W1,--no-as-needed 
all:
	rm -f egQ1.e
	make egQ1.e

egQ1.e: egQ1.cpp
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
asd.e: Ascend.cpp
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
scl.e: Scalar.cpp
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
trc.e: Trace.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd
