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

TEST(Tools, RDotR){

    Matrix A(4, 3);
    A.randomize();
    Matrix B(3, 4);
    B.randomize();
    ASSERT_EQ(A*B, RDotR(A, B));

}

TEST(Tools, RDotC){

    Matrix A(4, 3);
    A.randomize();
    Matrix B(CTYPE, 3, 4);
    B.randomize();
    ASSERT_EQ(A*B, RDotC(A, B));

}

TEST(Tools, CDotR){

    Matrix A(CTYPE, 4, 3);
    A.randomize();
    Matrix B(3, 4);
    B.randomize();
    ASSERT_EQ(A*B, CDotR(A, B));

}

TEST(Tools, CDotC){

    Matrix A(CTYPE, 4, 3);
    A.randomize();
    Matrix B(CTYPE, 3, 4);
    B.randomize();
    ASSERT_EQ(A*B, CDotC(A, B));

}

TEST(Tools, RAddR){
    
    Matrix A(10, 10);
    A.randomize();
    Matrix B(10, 10);
    B.orthoRand();
    double* elemA = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elemA, A.getElem(), sizeof(double)*A.elemNum());
    double* elemB = (double*)malloc(sizeof(double)*B.elemNum());
    memcpy(elemB, B.getElem(), sizeof(double)*B.elemNum());
    Matrix C0 = A;
    RAddR(C0, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0[i], elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(10, 10, true);
    Ad.randomize();
    double* elemAd = (double*)malloc(sizeof(double)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(), sizeof(double)*Ad.elemNum());
    Matrix C1 = Ad;
    RAddR(C1, B);
    int row = 0;
    for(size_t i = 0; i < C1.elemNum(); i++){
        if(i == 0 || i % (row*C1.col()+row) == 0){
            ASSERT_EQ(C1[i], elemAd[row] + elemB[i]);
            row++;
        }else
            ASSERT_EQ(C1[i], elemB[i]);
    }

    // If B is diagonal 
    Matrix Bd(10, 10, true);
    Bd.randomize();
    double* elemBd = (double*)malloc(sizeof(double)*Bd.elemNum());
    memcpy(elemBd, Bd.getElem(), sizeof(double)*Bd.elemNum());
    Matrix C2 = A;
    RAddR(C2, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2[i], elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2[i], elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = Ad;
    RAddR(C3, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3[i], elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);
}

TEST(Tools, CAddC){

    Matrix A(CTYPE, 10, 10);
    A.randomize();
    Matrix B(CTYPE, 10, 10);
    B.orthoRand();
    std::complex<double>* elemA = (std::complex<double>*)malloc(sizeof(std::complex<double>)*A.elemNum());
    memcpy(elemA, A.getElem(CTYPE), sizeof(std::complex<double>)*A.elemNum());
    std::complex<double>* elemB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*B.elemNum());
    memcpy(elemB, B.getElem(CTYPE), sizeof(std::complex<double>)*B.elemNum());
    Matrix C0 = A;
    CAddC(C0, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(CTYPE, 10, 10, true);
    Ad.randomize();
    std::complex<double>* elemAd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(CTYPE), sizeof(std::complex<double>)*Ad.elemNum());
    Matrix C1 = Ad;
    CAddC(C1, B);
    int row = 0;
    for(size_t i = 0; i < C1.elemNum(); i++){
        if(i == 0 || i % (row*C1.col()+row) == 0){
            ASSERT_EQ(C1(i), elemAd[row] + elemB[i]);
            row++;
        }else
            ASSERT_EQ(C1(i), elemB[i]);
    }

    // If B is diagonal 
    Matrix Bd(CTYPE, 10, 10, true);
    Bd.randomize();
    std::complex<double>* elemBd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Bd.elemNum());
    memcpy(elemBd, Bd.getElem(CTYPE), sizeof(std::complex<double>)*Bd.elemNum());
    Matrix C2 = A;
    CAddC(C2, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = Ad;
    CAddC(C3, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

TEST(Tools, RAddC){

    Matrix A(10, 10);
    A.randomize();
    Matrix B(CTYPE, 10, 10);
    B.orthoRand();
    double* elemA = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elemA, A.getElem(), sizeof(double)*A.elemNum());
    std::complex<double>* elemB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*B.elemNum());
    memcpy(elemB, B.getElem(CTYPE), sizeof(std::complex<double>)*B.elemNum());

    Matrix C0 = A;
    RAddC(C0, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(10, 10, true);
    Ad.randomize();
    double* elemAd = (double*)malloc(sizeof(double)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(), sizeof(double)*Ad.elemNum());
    Matrix C1 = Ad;
    RAddC(C1, B);
    int row = 0;
    for(size_t i = 0; i < C1.elemNum(); i++){
        if(i == 0 || i % (row*C1.col()+row) == 0){
            ASSERT_EQ(C1(i), elemAd[row] + elemB[i]);
            row++;
        }else
            ASSERT_EQ(C1(i), elemB[i]);
    }

    // If B is diagonal 
    Matrix Bd(CTYPE, 10, 10, true);
    Bd.randomize();
    std::complex<double>* elemBd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Bd.elemNum());
    memcpy(elemBd, Bd.getElem(CTYPE), sizeof(std::complex<double>)*Bd.elemNum());
    Matrix C2 = A;
    RAddC(C2, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = Ad;
    RAddC(C3, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

TEST(Tools, CAddR){

    Matrix A(CTYPE, 10, 10);
    A.randomize();
    Matrix B(10, 10);
    B.orthoRand();
    std::complex<double>* elemA = (std::complex<double>*)malloc(sizeof(std::complex<double>)*A.elemNum());
    memcpy(elemA, A.getElem(CTYPE), sizeof(std::complex<double>)*A.elemNum());
    double* elemB = (double*)malloc(sizeof(double)*B.elemNum());
    memcpy(elemB, B.getElem(), sizeof(double)*B.elemNum());

    Matrix C0 = A;
    CAddR(C0, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(CTYPE, 10, 10, true);
    Ad.randomize();
    std::complex<double>* elemAd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(CTYPE), sizeof(std::complex<double>)*Ad.elemNum());
    Matrix C1 = Ad;
    CAddR(C1, B);
    int row = 0;
    for(size_t i = 0; i < C1.elemNum(); i++){
        if(i == 0 || i % (row*C1.col()+row) == 0){
            ASSERT_EQ(C1(i), elemAd[row] + elemB[i]);
            row++;
        }else
            ASSERT_EQ(C1(i), elemB[i]);
    }

    // If B is diagonal 
    Matrix Bd(10, 10, true);
    Bd.randomize();
    double* elemBd = (double*)malloc(sizeof(double)*Bd.elemNum());
    memcpy(elemBd, Bd.getElem(), sizeof(double)*Bd.elemNum());
    Matrix C2 = A;
    CAddR(C2, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = Ad;
    CAddR(C3, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

