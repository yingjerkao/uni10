CXX       = icpc
CXXFLAGS += -O2 -std=c++11
INCFLAGS += -I /usr/local/uni10/include
LDFLAGS  += -L /usr/local/uni10/lib
LIBS     += -mkl -luni10


all: itebd_u1 genOP_u1

itebd_u1: itebd_u1.o func_itebd_u1.o func_evol_u1.o func_net_u1.o func_init_u1.o
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

genOP_u1: operator_u1.o
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $(LIBS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCFLAGS) $(LDFLAGS) $(LIBS) -c $^

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $^


.PHONY: clean init

clean:
	rm -f *.o 

init:
	rm -f *.o itebd_u1 genOP_u1
