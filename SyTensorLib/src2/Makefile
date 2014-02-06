OBJ := qnum.o bond.o block.o uni10_tools.o uni10_lapack.o matrix.o operate.o syTensor_elem.o syTensor.o network.o
LINK := -mkl=parallel
CC := icpc
UNI10_ROOT := uni10
SRC := .
BIN := bin

all: $(OBJ) libtensorF.a


doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

qnum.o: $(UNI10_ROOT)/datatype/lib/Qnum.cpp $(UNI10_ROOT)/datatype/Qnum.h
	$(CC) -c -I $(SRC) $< -o $@

bond.o: $(UNI10_ROOT)/data-structure/lib/Bond.cpp $(UNI10_ROOT)/data-structure/Bond.h
	$(CC) -c -I $(SRC) $< -o $@

block.o: $(UNI10_ROOT)/data-structure/lib/Block.cpp $(UNI10_ROOT)/data-structure/Block.h
	$(CC) -c -I $(SRC) $< -o $@

uni10_tools.o: $(UNI10_ROOT)/tools/lib/uni10_tools.cpp $(UNI10_ROOT)/tools/uni10_tools.h
	$(CC) -c -I $(SRC) $< -o $@

uni10_lapack.o: $(UNI10_ROOT)/numeric/lib/uni10_lapack.cpp $(UNI10_ROOT)/numeric/uni10_lapack.h
	$(CC) -c -I $(SRC) $< -o $@

matrix.o: $(UNI10_ROOT)/tensor-network/lib/Matrix.cpp $(UNI10_ROOT)/tensor-network/Matrix.h
	$(CC) -c -I $(SRC) $< -o $@

operate.o: $(UNI10_ROOT)/tensor-network/lib/operate.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c -I $(SRC) $< -o $@

syTensor_elem.o: $(UNI10_ROOT)/tensor-network/lib/SyTensor_elem.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c -I $(SRC) $< -o $@

syTensor.o: $(UNI10_ROOT)/tensor-network/lib/SyTensor.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c -I $(SRC) $< -o $@

network.o: $(UNI10_ROOT)/tensor-network/lib/Network.cpp $(UNI10_ROOT)/tensor-network/Network.h
	$(CC) -c -I $(SRC) $< -o $@

libtensorF.a: $(OBJ)
	rm -f ../lib/$@
	ar q ../lib/$@ $^

clean:
	rm -f *.o *.o *.e ../lib/*.a
