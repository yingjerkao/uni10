INC := ../../../SyTensorLib/src2/
MERA_INC := ../../MERALib/
LIB := ../../../SyTensorLib/src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

asd.e: Ascend.cpp
	icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)

flt.e: Filter.cpp
	icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)

main.e: main.cpp
	icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)


clean:
	rm -f *.o *.e asd
