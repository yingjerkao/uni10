OBJ := qnum.o bond.o block.o myLapack.o operate.o syTensor.o network.o matrix.o
LINK := -mkl=parallel

doc:
	doxygen Doxyfile
	cd ../docs/latex/; make pdf

all: $(OBJ)

qnum.o: Qnum.cpp Qnum.h
	icpc -c $< -o $@

bond.o: Bond.cpp Bond.h
	icpc -c $< -o $@

block.o: Block.cpp Block.h
	icpc -c $< -o $@

myLapack.o: myLapack.cpp myLapack.h
	icpc -c $< -mkl=parallel -o $@

operate.o: operate.cpp
	icpc -c $< -o $@

syTensor.o: SyTensor.cpp SyTensor.h
	icpc -c $< -o $@

matrix.o: Matrix.cpp Matrix.h
	icpc -c $< -o $@

network.o: Network.cpp Network.h
	icpc -c $< -o $@

clean:
	rm -f *.o *.e
