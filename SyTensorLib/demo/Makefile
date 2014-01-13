SRC := ../src/
INC := ../src/
LIB := ../lib/
OBJ := $(LIB)libtensorF.a
LINK := -mkl=parallel

mat.e: mat.cpp
	icpc -I$(INC) -L$(LIB) $(LINK) -o $@ $< $(OBJ)

syt_basic.e: sytensor_basic.cpp
	cd $(src); make all
	icpc -i$(inc) -l$(lib) $(link) -o $@ $< $(obj)


clean:
	rm -f *.o *.e
