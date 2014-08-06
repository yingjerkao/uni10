UNI10_ROOT := /Users/yunda/Research/TensorNetwork/uni10
INC := $(UNI10_ROOT)/include/
LIB := $(UNI10_ROOT)/lib/
CC:=g++
FLAGS:=-O3 -m64
#OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o $(UNI10_ROOT)/uni10_lapack.o $(UNI10_ROOT)/uni10_tools.o
#OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o

_pyUni10.so: pyUni10.i
	swig -c++ -python $<
	$(CC) -fPIC -I$(INC) -c pyUni10_wrap.cxx -I/usr/include/python2.4/
	$(CC) -I$(INC) -L$(LIB) $(FLAGS)  -shared  -o $@ -lblas -llapack -lm -luni10 pyUni10_wrap.o

uni10_mac: pyUni10.i
	swig -c++ -python $<
	$(CC) -c `python-config --cflags` -I$(INC) pyUni10_wrap.cxx
	$(CC) -I$(INC) -L$(LIB)  -bundle `python-config --ldflags` $(FLAGS) -o _pyUni10.so -lblas -llapack -lm -luni10 pyUni10_wrap.o
#g++ -shared $(LINK) $(OBJ) pyUni10_wrap.o -o $@
#g++ -shared $(OBJ) pyUni10_wrap.o -o $@ -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
#g++ $< -o $@ $(OBJ) -I$(INC) -I$(MERA_INC)  -DMKL_ILP64 -m64 -I$(MKLROOT)/include  -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -ldl -lpthread -lm



clean:
	rm -f *_wrap.* *.o *.so
