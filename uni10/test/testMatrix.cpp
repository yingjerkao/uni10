//
//  testBlock.cpp
//
//
//  Created by Ying-Jer Kao on 4/28/15.
//
//

#include <gtest/gtest.h>
#include <iostream>
#include <map>
#include "uni10.hpp"
#include <time.h>
#include <vector>
using namespace uni10;

TEST(Matrix, OperationAdd){

    // R + R
    Matrix A1(4, 3);
    A1.randomize();
    Matrix B1(4, 3);
    B1.orthoRand();
    Matrix Mc0 = A1;
    Mc0 += B1;
    Matrix Mc1 = A1 + B1;
    ASSERT_EQ(Mc0, A1+B1);
    for(size_t i = 0; i < A1.elemNum(); i++){
        ASSERT_EQ(Mc0[i], A1[i]+B1[i]);  
        ASSERT_EQ(Mc1[i], A1[i]+B1[i]);  
    }

    // C + R
    Matrix A2(CTYPE, 4, 3);
    A2.randomize();
    Matrix Md0 = A2;
    Md0 += B1;
    Matrix Md1 = B1 + A2;
    ASSERT_EQ(Md0, Md1);
    for(size_t i = 0; i < A2.elemNum(); i++){
        ASSERT_EQ(Md0(i), A2(i)+B1[i]);  
        ASSERT_EQ(Md1(i), A2(i)+B1[i]);  
    }

}
