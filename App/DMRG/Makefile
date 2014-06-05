INC := ../SyTensorLib/src2/
LIB := ../SyTensorLib/src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

mps3.e: MPS3.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
mps2.e: MPS2.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
mps.e: MPS.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
dmrg2.e: DMRG2.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
dmrg3.e: DMRG3.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
dmrg.e: DMRG.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd
