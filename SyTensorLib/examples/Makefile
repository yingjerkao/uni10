SRC := ../src/
INC := ../src/
LIB := ../src/
OBJ := $(LIB)qnum.o $(LIB)bond.o $(LIB)block.o $(LIB)myLapack.o $(LIB)operate.o $(LIB)syTensor.o $(LIB)network.o $(LIB)matrix.o
LINK := -mkl=parallel

mat.e: mat.cpp
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
net.e: net.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
test.e: test.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
asdnet.e: ascend_net.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
init.e: init_H.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
asd.e: ascend.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
rasd.e: rawAscend.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
dir.e: direct.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
demo.e: demo_s1.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
con.e: contract.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
run.e: main.cpp $(OBJ)
	cd $(SRC); make all
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)


clean:
	rm -f *.o *.e
