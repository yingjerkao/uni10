### BUILD
SRC_PATH=$(PWD)/../src
INCLUDE_PATH=$(PWD)/../src
TEST_SRC_PATH=$(PWD)/../test/src

BUILD_BIN_PATH=$(PWD)/bin
BUILD_TMP_PATH=$(PWD)/tmp
BUILD_LIB_PATH=$(PWD)/lib

### COMPILER 

CXX=g++
CXXFLAGS=-std=c++11 -DNDEBUG
INCLUDES=\
  -I$(INCLUDE_PATH)

OBJECTS=\
 $(BUILD_TMP_PATH)/Qnum.o \



### DO NOW
all: lib test

lib: $(OBJECTS)
	$(CXX) -shared -o $(BUILD_LIB_PATH)/libuni10.so $^ $(LIBRARIES)

test: $(BUILD_BIN_PATH)/Qnum.testCpp 


$(BUILD_TMP_PATH)/Qnum.o: $(SRC_PATH)/uni10/datatype/lib/Qnum.cpp $(SRC_PATH)/uni10/datatype/Qnum.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -fpic -c $<

$(BUILD_BIN_PATH)/Qnum.testCpp: $(TEST_SRC_PATH)/uni10/datatype/Qnum.testCpp.cpp $(lib)
	$(CXX) -I$(INCLUDE_PATH) -o $@ $< -L$(BUILD_LIB_PATH) -luni10

