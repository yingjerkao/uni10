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

TEST(Matrix, absMax){
    
    double elem1[12] = {9, 10, 29, -2, -100, 392, 33, -400, -12, -32, 12, 1};
    double elem2[12] = {9, 10, 29, -2, -100, 401, 33, -400, -12, -32, 12, 1};
    // R + R
    Matrix A1(4, 3, elem1);
    ASSERT_EQ(A1.absMax(), -400);
    Matrix A2(3, 4, elem2);
    ASSERT_EQ(A2.absMax(), 401);

}

TEST(Matrix, normalize){

    bool flag;
    Matrix A(3, 8);
    A.randomize();

    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), A.elemNum()*sizeof(double));
    double norm = A.norm(); 
    Matrix B = A.normalize();
    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / norm) < 1E-8;
        ASSERT_EQ(flag, true);
        flag = fabs(B.norm() - 1.) < 1E-8;
        ASSERT_EQ(flag, true);
    }

}

TEST(Matrix, absMaxNorm){

    bool flag;
    Matrix A(3, 8);
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

    Matrix B = A.absMaxNorm();

    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / absMax) < 1E-8;
        ASSERT_EQ(flag, true);
    }

}

TEST(Matrix, maxNorm){

    bool flag;
    Matrix A(3, 8);
    A.randomize();

    double* elem = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elem, A.getElem(), A.elemNum()*sizeof(double));
    
    double maxTmp = elem[0];
    for(size_t i = 0; i < A.elemNum(); i++)
        if(maxTmp < elem[i])
            maxTmp = elem[i];

    double max = A.max(); 
    ASSERT_EQ(maxTmp, max); 

    Matrix B = A.maxNorm();
    
    for(size_t i = 0; i < A.elemNum(); i++){
        flag = fabs(A[i] - B[i]) < 1E-8;
        ASSERT_EQ(flag, true); 
        flag = fabs(A[i] - elem[i] / max) < 1E-8;
        ASSERT_EQ(flag, true);
    }
}
