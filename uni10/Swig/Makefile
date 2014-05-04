UNI10_ROOT := ../../SyTensorLib/src2
OBJ := $(UNI10_ROOT)/qnum.o $(UNI10_ROOT)/bond.o

_pyUni10.so: pyUni10.i
	swig -c++ -python $<
	g++ -fPIC -I$(UNI10_ROOT) -c pyUni10_wrap.cxx -I/usr/include/python2.4/
	g++ -shared $(OBJ) pyUni10_wrap.o -o $@


clean:
	rm -f *_wrap.* *.o *.so
