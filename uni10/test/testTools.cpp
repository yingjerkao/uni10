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

TEST(Matrix, RDotR){

    Matrix A(4, 3);
    A.randomize();
    Matrix B(3, 4);
    B.randomize();
    ASSERT_EQ(A*B, RDotR(A, B));

}

TEST(Matrix, RDotC){

    Matrix A(4, 3);
    A.randomize();
    Matrix B(CTYPE, 3, 4);
    B.randomize();
    ASSERT_EQ(A*B, RDotC(A, B));

}

TEST(Matrix, CDotR){

    Matrix A(CTYPE, 4, 3);
    A.randomize();
    Matrix B(3, 4);
    B.randomize();
    ASSERT_EQ(A*B, CDotR(A, B));

}

TEST(Matrix, CDotC){

    Matrix A(CTYPE, 4, 3);
    A.randomize();
    Matrix B(CTYPE, 3, 4);
    B.randomize();
    ASSERT_EQ(A*B, CDotC(A, B));

}

TEST(Matrix, RAddR){
    
    Matrix A(10, 10);
    A.randomize();
    Matrix B(10, 10);
    B.orthoRand();
    double* elemA = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elemA, A.getElem(), sizeof(double)*A.elemNum());
    double* elemB = (double*)malloc(sizeof(double)*B.elemNum());
    memcpy(elemB, B.getElem(), sizeof(double)*B.elemNum());
    Matrix C0 = RAddR(A, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0[i], elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(10, 10, true);
    Ad.randomize();
    double* elemAd = (double*)malloc(sizeof(double)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(), sizeof(double)*Ad.elemNum());
    Matrix C1 = RAddR(Ad, B);
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
    Matrix C2 = RAddR(A, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2[i], elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2[i], elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = RAddR(Ad, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3[i], elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);
}

TEST(Matrix, CAddC){

    Matrix A(CTYPE, 10, 10);
    A.randomize();
    Matrix B(CTYPE, 10, 10);
    B.orthoRand();
    std::complex<double>* elemA = (std::complex<double>*)malloc(sizeof(std::complex<double>)*A.elemNum());
    memcpy(elemA, A.getElem(CTYPE), sizeof(std::complex<double>)*A.elemNum());
    std::complex<double>* elemB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*B.elemNum());
    memcpy(elemB, B.getElem(CTYPE), sizeof(std::complex<double>)*B.elemNum());

    Matrix C0 = CAddC(A, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(CTYPE, 10, 10, true);
    Ad.randomize();
    std::complex<double>* elemAd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(CTYPE), sizeof(std::complex<double>)*Ad.elemNum());
    Matrix C1 = CAddC(Ad, B);
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
    Matrix C2 = CAddC(A, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = CAddC(Ad, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

TEST(Matrix, RAddC){

    Matrix A(10, 10);
    A.randomize();
    Matrix B(CTYPE, 10, 10);
    B.orthoRand();
    double* elemA = (double*)malloc(sizeof(double)*A.elemNum());
    memcpy(elemA, A.getElem(), sizeof(double)*A.elemNum());
    std::complex<double>* elemB = (std::complex<double>*)malloc(sizeof(std::complex<double>)*B.elemNum());
    memcpy(elemB, B.getElem(CTYPE), sizeof(std::complex<double>)*B.elemNum());

    Matrix C0 = RAddC(A, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(10, 10, true);
    Ad.randomize();
    double* elemAd = (double*)malloc(sizeof(double)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(), sizeof(double)*Ad.elemNum());
    Matrix C1 = RAddC(Ad, B);
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
    Matrix C2 = RAddC(A, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = RAddC(Ad, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

TEST(Matrix, CAddR){

    Matrix A(CTYPE, 10, 10);
    A.randomize();
    Matrix B(10, 10);
    B.orthoRand();
    std::complex<double>* elemA = (std::complex<double>*)malloc(sizeof(std::complex<double>)*A.elemNum());
    memcpy(elemA, A.getElem(CTYPE), sizeof(std::complex<double>)*A.elemNum());
    double* elemB = (double*)malloc(sizeof(double)*B.elemNum());
    memcpy(elemB, B.getElem(), sizeof(double)*B.elemNum());

    Matrix C0 = CAddR(A, B);
    for(size_t i = 0; i < A.elemNum(); i++)
        ASSERT_EQ(C0(i), elemA[i] + elemB[i]);

    // If A is diagonal 
    Matrix Ad(CTYPE, 10, 10, true);
    Ad.randomize();
    std::complex<double>* elemAd = (std::complex<double>*)malloc(sizeof(std::complex<double>)*Ad.elemNum());
    memcpy(elemAd, Ad.getElem(CTYPE), sizeof(std::complex<double>)*Ad.elemNum());
    Matrix C1 = CAddR(Ad, B);
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
    Matrix C2 = CAddR(A, Bd);
    row = 0;
    for(size_t i = 0; i < C2.elemNum(); i++){
        if(i == 0 || i % (row*C2.col()+row) == 0){
            ASSERT_EQ(C2(i), elemA[i] + elemBd[row]);
            row++;
        }else
            ASSERT_EQ(C2(i), elemA[i]);
    }

    // If A adn B are diagonal 
    Matrix C3 = CAddR(Ad, Bd);
    ASSERT_EQ(C3.isDiag(), true);

    for(size_t i = 0; i < C3.elemNum(); i++)
        ASSERT_EQ(C3(i), elemAd[i] + elemBd[i]);

    free(elemA);
    free(elemB);
    free(elemAd);
    free(elemBd);

}

