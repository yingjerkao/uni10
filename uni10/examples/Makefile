SRC := ../src2/
INC := /opt/uni10/Library/uni10-0.1.0/include/
LIB := /opt/uni10/Library/uni10-0.1.0/lib/
MKL_LIB:=/opt/intel/mkl/lib/intel64
MKL_INC:=/opt/intel/mkl/include
#OBJ := $(LIB)libtensorF.a
#LINK := -mkl=parallel -luni10
OPT:=-W1,--no-as-needed 
all:
	rm -f egQ1.e
	make egQ1.e

egQ1.e: egQ1.cpp
	#cd $(SRC); make all
	#icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
	#g++ -I$(INC) -I$(MKL_INC) -fopenmp -L${LIB} -m64 -o $@ -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -luni10 $<
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
asd.e: Ascend.cpp
	#cd $(SRC); make all
	#icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
	#g++ -I$(INC) -I$(MKL_INC) -fopenmp -L${LIB} -m64 -o $@ -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -luni10 $<
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
scl.e: Scalar.cpp
	#cd $(SRC); make all
	#icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
	#g++ -I$(INC) -I$(MKL_INC) -fopenmp -L${LIB} -m64 -o $@ -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -luni10 $<
	g++ -I$(INC) -L${LIB} -m64 -o $@ -lblas -llapack -lm -luni10 $<
trc.e: Trace.cpp
	#cd $(SRC); make all
	icpc -I$(INC) $(LINK) -o $@ $< $(OBJ)
	#g++ $< -o $@ $(OBJ) $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a -I$(INC) -I$(MERA_INC) -L${MKL_PATH} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L$(MKL_LIB) -liomp5
	#g++ -I$(INC) -I$(MERA_INC) -I$(MKL_INC) -m64 -o $@ $< $(OBJ) $(MKL_LIB)/libmkl_intel_ilp64.a $(MKL_LIB)/libmkl_sequential.a $(MKL_LIB)/libmkl_core.a -lpthread -DMKL_ILP64 $(MKL_LIB1)/libirc.a $(MKL_LIB1)/libimf.a




clean:
	rm -f *.o *.e asd
