OBJ := qnum.o bond.o block.o myLapack.o operate.o syTensor.o
LINK := -mkl=parallel

demo: demo_s1.cpp $(OBJ)
	icpc -o $@ $< $(OBJ) $(LINK)
con: contract.cpp $(OBJ)
	icpc -o $@ $< $(OBJ) $(LINK)
run: main.cpp $(OBJ)
	icpc -o $@ $< $(OBJ) $(LINK)

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

clean:
	rm -f *.o
	rm -f run con demo
