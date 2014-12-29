#UNI10_ROOT := /home/yundasie/guni10
UNI10_ROOT := /Users/yunda/myLibrary/guni10
INC := ../src/
LIB := $(UNI10_ROOT)/lib
CC:=g++
NVCC:=nvcc
FLAGS:=-O3 -m64
#OBJ := $(LIB)/libuni10.a
CPU_OBJ := ../src/libuni10_cpu.a
GPU_OBJ := ../src/libuni10_gpu.a
CULA_INC := /usr/local/cula/include
CULA_LIB := /usr/local/cula/lib64
CULA_LINK := -lcula_core -lcula_lapack -lcudart -lcublas

unit.e: UniT.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
dmrg.e: fDMRG_U1.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
asd.e: Ascend.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
lan.e: Lanczos.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
gten.e: Ten.cpp
	$(NVCC) -I$(INC) -I$(CULA_INC) -L$(CULA_LIB) $(FLAGS) -o $@ -lblas -llapack -lm $(CULA_LINK) $< $(GPU_OBJ)
ten.e: Ten.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
gmul.e: MMul.cpp
	$(NVCC) -I$(INC) -I$(CULA_INC) -L$(CULA_LIB) $(FLAGS) -o $@ -lblas -llapack -lm $(CULA_LINK) $< $(GPU_OBJ)
mul.e: MMul.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)
gmat.e: Mat.cpp
	$(NVCC) -I$(INC) -I$(CULA_INC) -L$(CULA_LIB) $(FLAGS) -o $@ -lblas -llapack -lm $(CULA_LINK) $< $(GPU_OBJ)
mat.e: Mat.cpp
	$(CC) -I$(INC) $(FLAGS) -o $@ -lblas -llapack -lm  $< $(CPU_OBJ)

clean:
	rm -f *.o *.e asd
