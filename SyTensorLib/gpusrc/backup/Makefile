#OBJ := qnum.o bond.o block.o myTools.o myLapack.o operate.o syTensor_elem.o syTensor.o network.o matrix.o
OBJF_GPU := qnumF.o bond.o block.o myTools.o myLapack.o operate.o syTensor_elem.o syTensor.o network.o matrix.o
#OBJF_GPU := qnumF.o bond.o block.o myTools_g.o myLapack_g.o operate_g.o syTensor_elem_g.o syTensor_g.o network.o matrix_g.o
LIB := ../lib/
LINK := -mkl=parallel
CC := nvcc
NVCC := nvcc
MKL_LIB := /aracbox/intel/XE2011/sp1/mkl/lib/intel64
MKL_INC := /aracbox/intel/XE2011/sp1/mkl/include
NVCC_FLAG := -arch=sm_13

gpu: $(OBJF) libtensorF_GPU.a


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

myTools.o: myTools.cu myTools.h
	$(CC) -c $< -o $@ $(NVCC_FLAG)

myLapack.o: myLapack.cu myLapack.h
	$(CC) -c $< -o $@ $(NVCC_FLAG)

operate.o: operate.cpp SyTensor.h
	$(CC) -c $< -o $@

syTensor_elem.o: SyTensor_elem.cpp SyTensor.h myTools.h
	$(CC) -c $< -o $@

syTensor.o: SyTensor.cpp SyTensor.h
	$(CC) -c $< -o $@

matrix.o: Matrix.cpp Matrix.h
	$(CC) -c $< -o $@

network.o: Network.cpp Network.h
	$(CC) -c $< -o $@

libtensorF_GPU.a: $(OBJF_GPU)
	rm -f $(LIB)*.a
	ar q $(LIB)$@ $^

clean:
	rm -f *.o *.o *.e $(LIB)*.a
