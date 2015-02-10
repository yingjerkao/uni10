#####################
# Example Makefile
#####################
UNI10_ROOT:= /usr/local/uni10/ 
INC := $(UNI10_ROOT)/include/
LIB := $(UNI10_ROOT)/lib/
CC:=g++
FLAGS:=-O3 -m64
TARGETS:=egB1.e egB2.e egM1.e egM2.e egN1.e egQ1.e egQ2.e egU1.e egU2.e egU3.e
all: $(TARGETS)

$(TARGETS):%.e:%.cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ -lblas -llapack -lm -luni10 $<


.phony: clean
clean:
	rm -f *.o *.e 
