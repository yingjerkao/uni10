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

TEST(UniTensor,Constructor){

    int chi = 10; 
    Bond bdi_chi(BD_IN, chi); 
    UniTensor A;
    ASSERT_EQ(A.typeID(), 1);
    ASSERT_EQ(A.bondNum(), 0);

}


TEST(UniTensor,setElem){
    std::vector<Bond> bondsA(3, Bond(BD_OUT, 2));
    bondsA[0] = Bond(BD_IN, 3);
    UniTensor A(bondsA);
    A.randomize();
    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), sizeof(double)*A.elemNum());
    UniTensor B(bondsA);
    B.setRawElem(elem);
    UniTensor C(bondsA);
    C.setElem(elem);

    Matrix M(3, 4);
    M.setElem(elem);
    UniTensor D(bondsA);
    D.putBlock(M);
    for(size_t i = 0; i < A.elemNum(); i++){
        ASSERT_EQ(A.at(i), B.at(i));
        ASSERT_EQ(C.at(i), C.at(i));
        ASSERT_EQ(C.at(i), D.at(i));
    }
}

TEST(UniTensor,setRawElem){

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

TEST(UniTensor, normalize){

    std::vector<Bond> bonds(3, Bond(BD_OUT, 4));
    bonds[0] = Bond(BD_IN, 3);
    bool flag;
    UniTensor A(bonds);
    A.randomize();

    ASSERT_EQ(A.const_getBlock().norm(), A.norm());
    
    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), A.elemNum()*sizeof(double));
    double norm = A.norm(); 
    UniTensor B = A.normalize();
    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / norm) < 1E-8;
        ASSERT_EQ(flag, true);
        flag = fabs(B.norm() - 1.) < 1E-8;
        ASSERT_EQ(flag, true);
    }

}

TEST(UniTensor, absMaxNorm){

    std::vector<Bond> bonds(3, Bond(BD_OUT, 4));
    bonds[0] = Bond(BD_IN, 3);
    bool flag;
    UniTensor A(bonds);
    A.randomize();
    
    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), A.elemNum()*sizeof(double));

    double absMaxTmp = fabs(elem[0]);
    size_t idx = 0;
    for(size_t i = 0; i < A.elemNum(); i++)
        if(absMaxTmp < fabs(elem[i])){
            absMaxTmp = fabs(elem[i]) ;
            idx = i;
        }
    
    double absMax = A.absMax(); 
    ASSERT_EQ(absMaxTmp, absMax); 

    UniTensor B = A.absMaxNorm();

    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / absMax) < 1E-8;
        ASSERT_EQ(flag, true);
    }

}

TEST(UniTensor, maxNorm){

    std::vector<Bond> bonds(3, Bond(BD_OUT, 4));
    bonds[0] = Bond(BD_IN, 3);
    bool flag;
    UniTensor A(bonds);
    A.randomize();

    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), A.elemNum()*sizeof(double));
    
    double maxTmp = elem[0];
    for(size_t i = 0; i < A.elemNum(); i++)
        if(maxTmp < elem[i])
            maxTmp = elem[i];

    double max = A.max(); 
    ASSERT_EQ(maxTmp, max); 

    UniTensor B = A.maxNorm();
    
    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / max) < 1E-8;
        ASSERT_EQ(flag, true);
    }
}

