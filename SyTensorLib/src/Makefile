OBJ := qnum.o bond.o block.o myTools.o myLapack.o operate.o syTensor_elem.o syTensor.o network.o matrix.o
OBJF := qnumF.o bond.o block.o myTools.o myLapack.o matrix.o operate.o syTensor_elem.o syTensor.o network.o
OBJF_GPU := qnumF.o bond.o block.o myTools_g.o myLapack_g.o matrix.o operate.o syTensor_elem_g.o syTensor.o network.o
LINK := -mkl=parallel
CC := icpc
NVCC := nvcc
MKL_LIB := /aracbox/intel/XE2011/sp1/mkl/lib/intel64
MKL_INC := /aracbox/intel/XE2011/sp1/mkl/include

all: $(OBJF) libtensorF.a

gpu: $(OBJF_GPU) libtensorF_GPU.a


#all: $(OBJ) libtensor.a

doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

qnum.o: Qnum.cpp Qnum.h
	$(CC) -c $< -o $@

qnumF.o: QnumF.cpp QnumF.h
	$(CC) -c $< -o $@

bond.o: Bond.cpp Bond.h
	$(CC) -c $< -o $@

block.o: Block.cpp Block.h
	$(CC) -c $< -o $@

myTools.o: myTools.cpp myTools.h
	$(CC) -c $< -o $@

myLapack.o: myLapack.cpp myLapack.h
	$(CC) -c $< -o $@

matrix.o: Matrix.cpp Matrix.h
	$(CC) -c $< -o $@

operate.o: operate.cpp SyTensor.h
	$(CC) -c $< -o $@

syTensor_elem.o: SyTensor_elem.cpp SyTensor.h
	$(CC) -c $< -o $@

syTensor.o: SyTensor.cpp SyTensor.h
	$(CC) -c $< -o $@


network.o: Network.cpp Network.h
	$(CC) -c $< -o $@

myTools_g.o: myTools.cu myTools.h
	$(CC) -c $< -o $@

myLapack_g.o: myLapack.cu myLapack.h
	$(CC) -c $< -o $@

syTensor_elem_g.o: SyTensor_elem.cu SyTensor.h
	$(CC) -c $< -o $@

libtensorF.a: $(OBJF)
	rm -f ../lib/$@
	ar q ../lib/$@ $^
libtensorF_GPU.a: $(OBJF_GPU)
	rm -f ../lib/$@
	ar q ../lib/$@ $^

libtensor.a: $(OBJ)
	[ -d ../lib ] || mkdir -p ../lib
	ar q ../lib/$@ $^

clean:
	rm -f *.o *.o *.e ../lib/libtensorF.a
