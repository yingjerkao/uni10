<<<<<<< Updated upstream

SRCDIR = .
TENSORLIB = $(SRCDIR)/../../SyTensorLib/src/


INC = -I$(TENSORLIB) -I$(SRCDIR)
OBJ = $(TENSORLIB)qnumF.o $(TENSORLIB)bond.o $(TENSORLIB)block.o $(TENSORLIB)myLapack.o $(TENSORLIB)operate.o $(TENSORLIB)syTensor.o $(TENSORLIB)network.o $(TENSORLIB)matrix.o
LDFLAGS =


CXX = icpc
CXXFLAGS = -O2 -DNDEBUG $(INC)

LINKER = icpc
LIB = -lboost_program_options -mkl=parallel


itebd : itebd.o
	$(LINKER) $(LDFLAGS) $< -o $@ $(LIB) $(OBJ)

itebd.o : itebd.cpp itebd.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


itebd_info : itebd_info.o
	$(LINKER) $(LDFLAGS) $< -o $@ $(LIB) $(OBJ)

itebd_info.o : itebd_info.cpp itebd.h
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY : clean

clean :
	rm -f *.o
=======
UNI10_ROOT := /home/yundasie/uni10
#UNI10_ROOT := /Users/yunda/myLibrary/uni10
INC := $(UNI10_ROOT)/include
LIB := $(UNI10_ROOT)/lib
CC:=g++
FLAGS:=-O3 -m64
OBJ := $(LIB)/libuni10.a

fitebd.e: FiTEBD.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
fitebd_u1.e: FiTEBD_U1.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
itebd_u1f.e: iTEBD_U1_fixed.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
titebd_u1f.e: TiTEBD_U1_fixed.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<
itebd.e: iTEBD.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<

clean:
	rm -f *.o *.e asd
>>>>>>> Stashed changes
