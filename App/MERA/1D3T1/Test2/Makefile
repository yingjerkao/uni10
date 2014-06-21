#UNI10_ROOT := /home/yundasie/uni10
UNI10_ROOT := /Users/yunda/myLibrary/uni10
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
MERA_INC := ../../MERALib/
CC:=g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a

asd.e: Ascend2.cpp
	$(CC) -I$(INC) -I$(MERA_INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<

flt.e: Filter.cpp
	#cd $(SRC); make all
	icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)
	#g++ $< -o $@ $(OBJ) $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a -I$(INC) -I$(MERA_INC) -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
	#g++ -I$(INC) -I$(MERA_INC) -I$(MKL_INC) -m64 -o $@ $< $(OBJ) $(MKL_LIB)/libmkl_intel_ilp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -lpthread -DMKL_ILP64 $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a

asd_gpu.e: Ascend.cu
	#cd $(SRC); make gpu
	nvcc -o $@ $(OBJ_GPU) $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a $< -I$(INC) -I$(MERA_INC) -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5



clean:
	rm -f *.o *.e asd
