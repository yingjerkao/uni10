SRC := ../../SyTensorLib/src/
INC := ../../SyTensorLib/src/
LIB := ../../SyTensorLib/src/
MERAINC := ../MERALib
OBJ := $(LIB)qnumF.o $(LIB)bond.o $(LIB)block.o $(LIB)myLapack.o $(LIB)operate.o $(LIB)syTensor.o $(LIB)network.o $(LIB)matrix.o
LINK := -mkl=parallel

run: main.cpp
	cd $(SRC); make all
	icpc -I$(INC) -I$(MERAINC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)

asdnet.e: ascend_net.cpp
	icpc -o $@ $(LINK) -I$(TLIB) $< $(OBJ)

clean:
	rm -f *.o *.e run test
