INC := ../SyTensorLib/src2/
LIB := ../SyTensorLib/src2/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

tising_2d4s.e: TransIsing_2D4S.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
main.e: main.cpp
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)

clean:
	rm -f *.o *.e asd

