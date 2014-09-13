UNI10_ROOT := /home/yundasie/uni10
#UNI10_ROOT := /Users/yunda/myLibrary/uni10
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC := g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a

NVCC := nvcc
NVCCOPT := -arch=sm_13
CULA_INC := /usr/local/cula/include
CULA_LIB := /usr/local/cula/lib64
CULA_LINK := -lcula_core -lcula_lapack -lcudart -lcublas



OLD_UNI10_ROOT := ../../uni10/src
OLD_INC := $(OLD_UNI10_ROOT)
OLD_LIB := $(OLD_UNI10_ROOT)
OLD_OBJ := $(OLD_UNI10_ROOT)/libuni10_old.a

NEW_UNI10_ROOT := ../../guni10/src
NEW_INC := $(NEW_UNI10_ROOT)
NEW_LIB := $(NEW_UNI10_ROOT)
CPU_OBJ := $(NEW_UNI10_ROOT)/libuni10_cpu.a
GPU_OBJ := $(NEW_UNI10_ROOT)/libuni10_gpu.a


dmrg4.e: MPS4.cpp
	$(CC) -I$(NEW_INC) -L$(NEW_LIB) $(FLAGS) -o $@ -lblas -llapack -lm $< $(CPU_OBJ)
dmrg_gpu.e: DMRG_MPS.cpp
	$(NVCC) -I$(NEW_INC) -I$(CULA_INC) -L$(CULA_LIB) $(FLAGS) $(NVCCOPT) -o $@ -lblas -llapack -lm $(CULA_LINK) $< $(GPU_OBJ)

dmrg_cpu.e: DMRG_MPS.cpp
	$(CC) -I$(NEW_INC) -L$(NEW_LIB) $(FLAGS) -o $@ -lblas -llapack -lm $< $(CPU_OBJ)



clean:
	rm -f *.o *.e asd
