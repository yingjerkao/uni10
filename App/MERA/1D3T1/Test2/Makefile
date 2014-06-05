INC := ../../../SyTensorLib/src2/
MERA_INC := ../../MERALib/
LIB := ../../../SyTensorLib/src2/
OBJ := $(LIB)libtensorF.a
OBJ_GPU := $(LIB)libtensorF_GPU.a
LINK := -mkl=parallel
MKL_PATH := /aracbox/intel/XE2011/sp1/mkl
MKL_LIB := /aracbox/intel/XE2011/sp1/mkl/lib/intel64
MKL_INC := /aracbox/intel/XE2011/sp1/mkl/include
MKL_LIB1 := /aracbox/intel/XE2011/sp1/composer_xe_2011_sp1.7.256/compiler/lib/intel64

asd.e: Ascend2.cpp
	#cd $(SRC); make all
	#icpc -I$(INC) -I$(MERA_INC) $(LINK) -o $@ $< $(OBJ)
	#g++ $< -o $@ $(OBJ) -I$(INC) -I$(MERA_INC) -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
	g++ $< -o $@ $(OBJ) $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a -I$(INC) -I$(MERA_INC) -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
	#g++ -I$(INC) -I$(MERA_INC) -I$(MKL_INC) -m64 -o $@ $< $(OBJ) $(MKL_LIB)/libmkl_intel_ilp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -lpthread -DMKL_ILP64 $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a
	#g++ $< -o $@ $(OBJ) -I$(INC) -I$(MERA_INC)  -DMKL_ILP64 -m64 -I$(MKLROOT)/include  -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl -lpthread -lm

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
