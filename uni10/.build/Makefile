### BUILD
SRC_PATH=$(PWD)/../src
INCLUDE_PATH=$(PWD)/../src

BUILD_TMP_PATH=$(PWD)/tmp
BUILD_LIB_PATH=$(PWD)/lib

### COMPILER 

CXX=g++
CXXFLAGS=-std=c++11
INCLUDES=\
  -I$(INCLUDE_PATH)

OBJECTS=\
 $(BUILD_TMP_PATH)/Qnum.o

### DO NOW
all: lib

lib: $(OBJECTS)
	$(CXX) -shared -o $(BUILD_LIB_PATH)/libuni10.so $^ $(LIBRARIES)

$(BUILD_TMP_PATH)/Qnum.o: $(SRC_PATH)/uni10/datatype/lib/Qnum.cpp $(SRC_PATH)/uni10/datatype/Qnum.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -fpic -c $<
