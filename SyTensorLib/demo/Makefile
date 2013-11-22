SRC := ../src/
INC := ../src/
LIB := ../src/
OBJ := $(LIB)qnum.o $(LIB)bond.o $(LIB)block.o $(LIB)myLapack.o $(LIB)operate.o $(LIB)syTensor.o $(LIB)network.o $(LIB)matrix.o
LINK := -mkl=parallel

syt_basic.e: SyTensor_basic.cpp
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)


clean:
	rm -f *.o *.e
