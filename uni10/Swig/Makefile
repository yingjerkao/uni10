UNI10_ROOT := ../../SyTensorLib/src2
LINK := -mkl=parallel
MKL_PATH := /aracbox/intel/XE2011/sp1/mkl
MKL_LIB := /aracbox/intel/XE2011/sp1/mkl/lib/intel64
MKL_LIB1 := /aracbox/intel/XE2011/sp1/composer_xe_2011_sp1.7.256/compiler/lib/intel64
OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o $(UNI10_ROOT)/matrix.o $(UNI10_ROOT)/uni10_tools.o $(UNI10_ROOT)/uni10_lapack.o
#OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o $(UNI10_ROOT)/uni10_lapack.o $(UNI10_ROOT)/uni10_tools.o
#OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o 

_pyUni10.so: pyUni10.i
	swig -c++ -python $<
	g++ -fPIC -I$(UNI10_ROOT) -c pyUni10_wrap.cxx -I/usr/include/python2.4/
	#g++ -shared $(OBJ) pyUni10_wrap.o -o $@
	g++ -shared -o $@ $(OBJ) pyUni10_wrap.o $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
	

#g++ -shared $(LINK) $(OBJ) pyUni10_wrap.o -o $@
#g++ -shared $(OBJ) pyUni10_wrap.o -o $@ -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
#g++ $< -o $@ $(OBJ) -I$(INC) -I$(MERA_INC)  -DMKL_ILP64 -m64 -I$(MKLROOT)/include  -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl -lpthread -lm



clean:
	rm -f *_wrap.* *.o *.so
