INC := /opt/uni10/Library/uni10-0.1.0/include
LIB := /opt/uni10/Library/uni10-0.1.0/lib
#UNI10:=libtensorF.a
LINK := -mkl=parallel -luni10
CXX=icpc

%.o: %.cpp
	$(CXX) -I$(INC) -c -g $^  
asd.e: ascend.o 
	$(CXX) -I$(INC) -L$(LIB) $(LINK) -g -o $@ $< $(UNI10)
.PHONY: clean 
clean:
	rm -f *.o *.e
