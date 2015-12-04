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

/*
 class Block{
	public:
 Block();
 Block(const Block& _b);
 Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
 virtual ~Block();


 double operator[](size_t idx)const;
 double at(size_t i, size_t j)const;
 double* getElem()const;
 Matrix getDiag()const;
 void save(const std::string& fname)const;
 std::vector<CMatrix> eig()const;
 std::vector<Matrix> eigh()const;
 std::vector<Matrix> svd()const;
 size_t lanczosEigh(double& E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15)const;

 double sum()const;
 friend Matrix operator*(const Block& Ma, const Block& Mb);
 friend Matrix operator*(double a, const Block& Ma);
 friend Matrix operator*(const Block& Ma, double a);
 friend CMatrix operator*(const std::complex<double>& a, const Block& Ma);
 friend CMatrix operator*(const Block& Ma, const std::complex<double>& a);
 friend Matrix operator+(const Block& Ma, const Block& Mb);
 friend bool operator==(const Block& m1, const Block& m2);
 friend class UniTensor;
 friend class CUniTensor;
 friend class CBlock;
 friend class Matrix;
 friend class CMatrix;
 friend std::ostream& operator<< (std::ostream& os, const Block& b);
 //friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
	protected:
 double* m_elem;
 size_t Rnum;		//number of rows of the block
 size_t Cnum;		//number of columns of the block
 bool diag;
 bool ongpu;
 };

*/


TEST(Block, DefaultConstructor){

    Block blk;
    EXPECT_EQ(blk.row(),0);
    EXPECT_EQ(blk.col(),0);
    EXPECT_FALSE(blk.isDiag());
    EXPECT_EQ(blk.elemNum(),0);
}

TEST(Block,Constructor){

    Block blk(3,5);

    EXPECT_EQ(blk.row(),3);
    EXPECT_EQ(blk.col(),5);
    EXPECT_FALSE(blk.isDiag());
    EXPECT_EQ(blk.elemNum(),15);

}

TEST(Block,CopyConstructor){

    Block blk1(5,7);
    Block blk2(blk1);

    EXPECT_EQ(5,blk2.row());
    EXPECT_EQ(7,blk2.col());

    EXPECT_EQ(35, blk2.elemNum());

    EXPECT_FALSE(blk2.isDiag());

}

TEST(Block,ConstructorDia){

    Block blk(3,5,true);

    EXPECT_EQ(3,blk.row());
    EXPECT_EQ(5,blk.col());
    EXPECT_TRUE(blk.isDiag());
    EXPECT_EQ(3,blk.elemNum());

}

TEST(Block,CopyConstructorDia){

    Block blk1(5,7,true);
    Block blk2(blk1);

    EXPECT_EQ(5,blk2.row());
    EXPECT_EQ(7,blk2.col());
    EXPECT_TRUE(blk2.isDiag());
    EXPECT_EQ(5,blk2.elemNum());

}


TEST(Block,getElem){
    /*Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    UniTensor T(bdlist,"T");
*/

    Matrix T(5,5);
    T.randomize();

    Block blk1=T;
    //Block blk1=T.const_getBlock();

    EXPECT_TRUE(blk1.getElem()!=NULL);


}

TEST(Block,getDiag){
    Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    double heisenberg_s1[] = \
    {1, 0, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 0, 1, 0, 0, 0, 0, 0,\
        0, 0,-1, 0, 1, 0, 0, 0, 0,\
        0, 1, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 1, 0, 0, 0, 1, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 1, 0,\
        0, 0, 0, 0, 1, 0,-1, 0, 0,\
        0, 0, 0, 0, 0, 1, 0, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 0, 1\
    };
    // The block does not allocate memory, and the only way is through UniTensor::const_getBlock()

    UniTensor H(bdlist,"H");
    H.setRawElem(heisenberg_s1);
    //std::cout<< H ;

    Block blk1,blk2;

    blk1=H.const_getBlock();

    //std::cout<< blk1;
    blk2=blk1.getDiag();

    //std::cout << blk2;
    EXPECT_DOUBLE_EQ(-1.0,blk2[0]);
    EXPECT_DOUBLE_EQ(0.0,blk2[1]);
    EXPECT_DOUBLE_EQ(-1.0,blk2[2]);

    EXPECT_DOUBLE_EQ(-2.0,blk2.trace());

}


TEST(Block,inverse){
    Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    // The block does not allocate memory, and the only way is through
    // UniTensor::const_getBlock()
    double heisenberg_s1[] = \
    {1, 0, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 0, 1, 0, 0, 0, 0, 0,\
        0, 0,-1, 0, 1, 0, 0, 0, 0,\
        0, 1, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 1, 0, 0, 0, 1, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 1, 0,\
        0, 0, 0, 0, 1, 0,-1, 0, 0,\
        0, 0, 0, 0, 0, 1, 0, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 0, 1\
    };
    UniTensor H(bdlist,"H");
    H.orthoRand();
    Block blk1=H.const_getBlock();
    Matrix I,blk2,blk3;
    blk2=blk1.inverse();
    I.resize(blk2.row(),blk2.col()).identity();
    blk3=blk2*blk1;
    blk3+=(-1.0)*I;
    EXPECT_NEAR(0.0,blk3.max(), 1e-10);
}

TEST(Block, trace){
    Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    // The block does not allocate memory, and the only way is through
    // UniTensor::const_getBlock()
    double heisenberg_s1[] = \
    {1, 0, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 0, 1, 0, 0, 0, 0, 0,\
        0, 0,-1, 0, 1, 0, 0, 0, 0,\
        0, 1, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 1, 0, 0, 0, 1, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 1, 0,\
        0, 0, 0, 0, 1, 0,-1, 0, 0,\
        0, 0, 0, 0, 0, 1, 0, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 0, 1\
    };
    UniTensor H(bdlist,"H");
    H.setRawElem(heisenberg_s1);
    Block blk1=H.const_getBlock();
    EXPECT_DOUBLE_EQ(-2.0,blk1.trace());
}

TEST(Block, norm){
    Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    // The block does not allocate memory, and the only way is through
    // UniTensor::const_getBlock()
    double heisenberg_s1[] = \
    {1, 0, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 0, 1, 0, 0, 0, 0, 0,\
        0, 0,-1, 0, 1, 0, 0, 0, 0,\
        0, 1, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 1, 0, 0, 0, 1, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 1, 0,\
        0, 0, 0, 0, 1, 0,-1, 0, 0,\
        0, 0, 0, 0, 0, 1, 0, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 0, 1\
    };
    UniTensor H(bdlist,"H");
    H.setRawElem(heisenberg_s1);
    Block blk1=H.const_getBlock();
    EXPECT_DOUBLE_EQ(sqrt(6.0),blk1.norm());
}

TEST(Block, sum){
    Qnum q0(0),q1(1),q_1(-1);
    std::vector<Qnum> qnums={q1,q0,q_1};
    Bond bd_out(BD_OUT, qnums),bd_in(BD_IN,qnums);
    std::vector<Bond> bdlist={bd_in,bd_in,bd_out,bd_out};
    // The block does not allocate memory, and the only way is through
    // UniTensor::const_getBlock()
    double heisenberg_s1[] = \
    {1, 0, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 0, 1, 0, 0, 0, 0, 0,\
        0, 0,-1, 0, 1, 0, 0, 0, 0,\
        0, 1, 0, 0, 0, 0, 0, 0, 0,\
        0, 0, 1, 0, 0, 0, 1, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 1, 0,\
        0, 0, 0, 0, 1, 0,-1, 0, 0,\
        0, 0, 0, 0, 0, 1, 0, 0, 0,\
        0, 0, 0, 0, 0, 0, 0, 0, 1\
    };
    UniTensor H(bdlist,"H");
    H.setRawElem(heisenberg_s1);
    Block blk1=H.const_getBlock();
    EXPECT_DOUBLE_EQ(2.0,blk1.sum());
}


