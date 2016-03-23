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

TEST(UniTensor,DefaultConstructor){

    UniTensor A;
    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(A.bondNum(), 0);

}

TEST(UniTensor, putBlock){
    
    Matrix MA(4, 9);  
    MA.randomize();
    ASSERT_EQ(MA.typeID(), 1);
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);
    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(A.getBlock(), MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A[i], MA[i]); 

    MA.assign(CTYPE, 4, 9);
    MA.randomize();
    ASSERT_EQ(MA.typeID(), 2);
    A.putBlock(MA);
    ASSERT_EQ(A.typeID(), 2);
    ASSERT_EQ(A.getBlock(), MA);
    
    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A(i), MA(i)); 

}

TEST(UniTensor, at_bidx){
    
    Matrix MA(4, 9);  
    MA.randomize();
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A.at(i), MA[i]);

}

TEST(UniTensor, at_rflag_idx){

    Matrix MA(4, 9);  
    MA.randomize();
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A.at(RTYPE, i), MA[i]);

}

TEST(UniTensor, at_cflag_idx){

    Matrix MA(CTYPE, 4, 9);  
    MA.randomize();
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A.at(CTYPE, i), MA(i));

}

TEST(UniTensor, Brackets){

    Matrix MA(4, 9);  
    MA.randomize();
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A[i], MA[i]);


}

TEST(UniTensor, Parentheses){

    Matrix MA(CTYPE, 4, 9);  
    MA.randomize();
    std::vector<Bond> bondsA(4, Bond(BD_OUT, 3));
    bondsA[0] = Bond(BD_IN, 2);
    bondsA[1] = bondsA[0];
    UniTensor A(bondsA);
    A.putBlock(MA);

    for(size_t i = 0; i < MA.elemNum(); i++)
        ASSERT_EQ(A(i), MA(i));

}

TEST(UniTensor,Constructor){

    int chi = 10; 
    Bond bdi_chi(BD_IN, chi); 
    UniTensor A;
    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(A.bondNum(), 0);

}


