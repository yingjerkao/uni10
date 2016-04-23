//
//  Test_UniTensor.cpp
//  
//
//  Created by Ying-Jer Kao on 3/20/15.
//
//
#include <gtest/gtest.h>
#include <iostream>
#include <map>
#include "uni10.hpp"
#include <time.h>
#include <vector>
using namespace uni10;


TEST(Network,SimpleContract){
    
    Matrix MA(3, 4);
    MA.randomize();
    Matrix MB(4, 3);
    MB.randomize();
    std::vector<Bond> bondsA(3, Bond(BD_OUT, 2));
    bondsA[0] = Bond(BD_IN, 3);
    UniTensor A(bondsA);
    UniTensor B = A;
    A.putBlock(MA);
    B.transpose();
    B.putBlock(MB);

    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(B.typeID(), 1);

    Network Simple_net("./Simple.net");
    Simple_net.putTensor(0, A);
    Simple_net.putTensor(1, B);
    UniTensor C = Simple_net.launch();

    ASSERT_EQ(C.getBlock(), MA*MB);
    ASSERT_EQ(C.typeID(), 1);

    MA.assign(CTYPE, 3, 4);
    MA.randomize();
    MB.assign(CTYPE, 4, 3);
    MB.randomize();
    A.putBlock(MA);
    B.putBlock(MB);

    ASSERT_EQ(A.typeID(), 2);
    ASSERT_EQ(B.typeID(), 2);

    Simple_net.putTensor(0, A);
    Simple_net.putTensor(1, B);
    C = Simple_net.launch();

    ASSERT_EQ(C.getBlock(), MA*MB);
    ASSERT_EQ(C.typeID(), 2);

    MA.assign(RTYPE, 3, 4);
    MA.randomize();
    MB.assign(CTYPE, 4, 3);
    MB.randomize();
    A.putBlock(MA);
    B.putBlock(MB);

    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(B.typeID(), 2);

    Simple_net.putTensor(0, A);
    Simple_net.putTensor(1, B);
    C = Simple_net.launch();

    ASSERT_EQ(C.getBlock(), MA*MB);
    ASSERT_EQ(C.typeID(), 2);

    MA.assign(CTYPE, 3, 4);
    MA.randomize();
    MB.assign(RTYPE, 4, 3);
    MB.randomize();
    A.putBlock(MA);
    B.putBlock(MB);

    ASSERT_EQ(A.typeID(), 2);
    ASSERT_EQ(B.typeID(), 1);

    Simple_net.putTensor(0, A);
    Simple_net.putTensor(1, B);
    C = Simple_net.launch();

    ASSERT_EQ(C.getBlock(), MA*MB);
    ASSERT_EQ(C.typeID(), 2);

}
