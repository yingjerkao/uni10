SRC := ../src2/
INC := ../src2/
LIB := ../src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

all:
	rm -f egM3.e
	make egM3.e

egN1.e: egN1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egM2.e: egM2.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egM1.e: egM1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egM3.e: egM3.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egU1.e: egU1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egU2.e: egU2.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egU3.e: egU3.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egB1.e: egB1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egQ1.e: egQ1.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
egQ2.e: egQ2.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
scl.e: Scalar.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
trc.e: Trace.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)




clean:
	rm -f *.o *.e asd
