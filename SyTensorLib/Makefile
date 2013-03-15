run: main.cpp Qnum.h Qnum.cpp Bond.h Bond.cpp Block.h SyTensor.h SyTensor.cpp
	g++ -c Qnum.cpp
	g++ -c Block.cpp
	g++ -c Bond.cpp
	g++ -c SyTensor.cpp
	g++ -o $@ $< Qnum.o Bond.o SyTensor.o Block.o
