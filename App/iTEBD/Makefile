
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
