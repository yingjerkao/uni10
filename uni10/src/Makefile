OBJ := qnum.o bond.o block.o myLapack.o operate.o syTensor.o network.o matrix.o
OBJF := qnumF.o bond.o block.o myLapack.o operate.o syTensor.o network.o matrix.o
LINK := -mkl=parallel
CC := icpc

all: $(OBJF) libtensorF.a

#all: $(OBJ) libtensor.a

doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

qnum.o: Qnum.cpp Qnum.h
	$(CC) -c $< -o $@

qnumF.o: QnumF.cpp QnumF.h
	$(CC) -c $< -o $@

bond.o: Bond.cpp Bond.h
	$(CC) -c $< -o $@

block.o: Block.cpp Block.h
	$(CC) -c $< -o $@

myLapack.o: myLapack.cpp myLapack.h
	$(CC) -c $< -mkl=parallel -o $@

operate.o: operate.cpp
	$(CC) -c $< -o $@

syTensor.o: SyTensor.cpp SyTensor.h
	$(CC) -c $< -o $@

matrix.o: Matrix.cpp Matrix.h
	$(CC) -c $< -o $@

network.o: Network.cpp Network.h
	$(CC) -c $< -o $@

libtensorF.a: $(OBJF)
	ar q ../lib/$@ $^

libtensor.a: $(OBJ)
	[ -d ../lib ] || mkdir -p ../lib
	ar q ../lib/$@ $^

clean:
	rm -f *.o *.e
