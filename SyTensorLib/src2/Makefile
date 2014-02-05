OBJ := qnumF.o bond.o block.o myTools.o myLapack.o matrix.o operate.o syTensor_elem.o syTensor.o network.o
LINK := -mkl=parallel
CC := icpc
UNI10_ROOT := uni10
BIN := bin

all: $(OBJ) libtensorF.a


doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

qnumF.o: $(UNI10_ROOT)/datatype/lib/QnumF.cpp $(UNI10_ROOT)/datatype/QnumF.h
	$(CC) -c $< -o $@

bond.o: $(UNI10_ROOT)/datatype/lib/Bond.cpp $(UNI10_ROOT)/datatype/Bond.h
	$(CC) -c $< -o $@

block.o: $(UNI10_ROOT)/data-structure/lib/Block.cpp $(UNI10_ROOT)/data-structure/Block.h
	$(CC) -c $< -o $@

myTools.o: $(UNI10_ROOT)/tensor-network/lib/myTools.cpp $(UNI10_ROOT)/tensor-network/myTools.h
	$(CC) -c $< -o $@

myLapack.o: $(UNI10_ROOT)/numeric/lib/myLapack.cpp $(UNI10_ROOT)/numeric/myLapack.h
	$(CC) -c $< -o $@

matrix.o: $(UNI10_ROOT)/tensor-network/lib/Matrix.cpp $(UNI10_ROOT)/tensor-network/Matrix.h
	$(CC) -c $< -o $@

operate.o: $(UNI10_ROOT)/tensor-network/lib/operate.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c $< -o $@

syTensor_elem.o: $(UNI10_ROOT)/tensor-network/lib/SyTensor_elem.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c $< -o $@

syTensor.o: $(UNI10_ROOT)/tensor-network/lib/SyTensor.cpp $(UNI10_ROOT)/tensor-network/SyTensor.h
	$(CC) -c $< -o $@


network.o: $(UNI10_ROOT)/tensor-network/lib/Network.cpp $(UNI10_ROOT)/tensor-network/Network.h
	$(CC) -c $< -o $@

libtensorF.a: $(OBJ)
	rm -f ../lib/$@
	ar q ../lib/$@ $^

clean:
	rm -f *.o *.o *.e ../lib/*.a
