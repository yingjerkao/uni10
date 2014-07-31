CPU_OBJ := qnum.o bond.o block.o uni10_tools.o uni10_tools_cpu.o uni10_lapack_cpu.o matrix.o operate.o UniTensor_elem.o UniTensor.o network.o vecadd.o
GPU_OBJ := qnum.o bond.o block.o uni10_tools.o uni10_tools_gpu.o uni10_dgemm.o uni10_lapack_gpu.o matrix.o operate.o UniTensor_elem.o UniTensor.o network.o vecadd.o
LINK := -mkl=parallel
CC := g++
UNI10_ROOT := uni10
SRC := .
BIN := bin
OPT := -std=c++0x -O3
NVCC := nvcc
NVCCOPT := -arch=sm_13
CULA_INC := /usr/local/cula/include
CULA_LIB := /usr/local/cula/lib64

gpu: $(GPU_OBJ) libuni10_gpu.a
cpu: $(CPU_OBJ) libuni10_cpu.a


doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

qnum.o: $(UNI10_ROOT)/datatype/lib/Qnum.cpp $(UNI10_ROOT)/datatype/Qnum.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

bond.o: $(UNI10_ROOT)/data-structure/lib/Bond.cpp $(UNI10_ROOT)/data-structure/Bond.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

block.o: $(UNI10_ROOT)/data-structure/lib/Block.cpp $(UNI10_ROOT)/data-structure/Block.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

uni10_tools.o: $(UNI10_ROOT)/tools/lib/uni10_tools.cpp $(UNI10_ROOT)/tools/uni10_tools.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

uni10_tools_cpu.o: $(UNI10_ROOT)/tools/lib/uni10_tools_cpu.cpp $(UNI10_ROOT)/tools/uni10_tools.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

uni10_tools_gpu.o: $(UNI10_ROOT)/tools/lib/uni10_tools_gpu.cu $(UNI10_ROOT)/tools/uni10_tools.h
	$(NVCC) -c $(NVCCOPT) -I $(SRC) $< -o $@

uni10_dgemm.o: $(UNI10_ROOT)/numeric/lib/uni10_dgemm.cu $(UNI10_ROOT)/numeric/uni10_lapack.h
	$(NVCC) -c $(NVCCOPT) -I $(SRC) -I$(CULA_INC) $< -o $@

uni10_lapack_cpu.o: $(UNI10_ROOT)/numeric/lib/uni10_lapack_cpu.cpp $(UNI10_ROOT)/numeric/uni10_lapack.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

uni10_lapack_gpu.o: $(UNI10_ROOT)/numeric/lib/uni10_lapack_gpu.cu $(UNI10_ROOT)/numeric/uni10_lapack.h
	$(NVCC) -c $(NVCCOPT) -I $(SRC) -I$(CULA_INC) $< -o $@

matrix.o: $(UNI10_ROOT)/tensor-network/lib/Matrix.cpp $(UNI10_ROOT)/tensor-network/Matrix.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

operate.o: $(UNI10_ROOT)/tensor-network/lib/operate.cpp $(UNI10_ROOT)/tensor-network/UniTensor.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

UniTensor_elem.o: $(UNI10_ROOT)/tensor-network/lib/UniTensor_elem.cpp $(UNI10_ROOT)/tensor-network/UniTensor.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

UniTensor.o: $(UNI10_ROOT)/tensor-network/lib/UniTensor.cpp $(UNI10_ROOT)/tensor-network/UniTensor.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

network.o: $(UNI10_ROOT)/tensor-network/lib/Network.cpp $(UNI10_ROOT)/tensor-network/Network.h
	$(CC) -c $(OPT) -I $(SRC) $< -o $@

vecadd.o: $(UNI10_ROOT)/tensor-network/lib/vecadd.cu $(UNI10_ROOT)/tensor-network/Network.h
	$(NVCC) -c $(NVCCOPT) -I $(SRC) $< -o $@

libuni10_cpu.a: $(CPU_OBJ)
	ar rcs ./$@ $^

libuni10_gpu.a: $(GPU_OBJ)
	ar rcs ./$@ $^

clean:
	rm -f *.o *.o *.e ./*.a
