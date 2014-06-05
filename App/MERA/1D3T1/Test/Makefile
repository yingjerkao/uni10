SRC := ../../../SyTensorLib/src/
INC := ../../../SyTensorLib/src/
MERA_INC := ../../MERALib/
LIB := ../../../SyTensorLib/lib/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel



asds.e: AscendS.cpp
	icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)
dess.e: DescendS.cpp
	cd $(SRC); make all
	icpc -I$(INC) -I$(MERA_INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)
asds1: AscendS1.cpp
	cd $(SRC); make all
	icpc -I$(INC) -I$(MERA_INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)

asd: Ascend.cpp
	cd $(SRC); make all
	icpc -I$(INC) -I$(MERA_INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd
