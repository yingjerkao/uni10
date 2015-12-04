//
//  Test_Qnum.cpp
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
TEST(Qnum,DefaultConstructor){
    // Qnum(int _U1=0, parityType _prt=PRT_EVEN);
    
    Qnum q;
    EXPECT_EQ(0,q.U1());
    EXPECT_EQ(PRT_EVEN, q.prt());
    EXPECT_FALSE(q.isFermionic());
}



TEST(Qnum,DefaultFermionConstructor){
    // Qnum(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);
    Qnum q(PRTF_ODD);
    EXPECT_EQ(PRTF_ODD, q.prtF());
    EXPECT_EQ(0,q.U1());
    EXPECT_EQ(PRT_EVEN,q.prt());
    EXPECT_TRUE(q.isFermionic());
}

TEST(Qnum,Constructor){
    EXPECT_THROW({
        Qnum q1(15000);
    }, std::exception);
}


TEST(Qnum,CopyConstructor){
    // Qnum(const Qnum& _q);
    Qnum q1(1);
    
    Qnum q2(q1);
    
    EXPECT_EQ(q1,q2);
    
}
TEST(Qnum,isFermionic){
    // Qnum(int _U1=0, parityType _prt=PRT_EVEN);
    Qnum q;
  
    EXPECT_EQ(0,q.U1());
    EXPECT_EQ(PRT_EVEN,q.prt());
// PRTF_ODD is defined, so fermionic is set true
    EXPECT_TRUE(q.isFermionic()); 
}

TEST(Qnum,assign){
    Qnum q;
    q.assign(3);
    EXPECT_EQ(3,q.U1());
    EXPECT_EQ(PRT_EVEN,q.prt());
    
    q.assign(0,PRT_ODD);
    EXPECT_EQ(PRT_ODD,q.prt());
    EXPECT_EQ(0,q.U1());
    
    q.assign(-1);
    EXPECT_EQ(PRT_EVEN,q.prt()); // Default value PRT_EVEN assigned
    EXPECT_EQ(-1,q.U1());
    
    Qnum qf;
    
    qf.assign(PRTF_EVEN, 1, PRT_ODD);
    
    EXPECT_EQ(qf.prtF(),PRTF_EVEN);
    EXPECT_EQ(qf.U1(),1);
    EXPECT_EQ(qf.prt(), PRT_ODD);

    qf.assign(PRTF_ODD);
    
    EXPECT_EQ(qf.prtF(), PRTF_ODD);
    EXPECT_EQ(qf.prt(),PRT_EVEN); // Default values are also implicitly assigned
    EXPECT_EQ(qf.U1(),0);

}

TEST(Qnum,hash){
    Qnum q(PRTF_ODD,100,PRT_EVEN);
    EXPECT_EQ(q.hash(), 401);
}


TEST(Qnum, Operation_Negation){
    Qnum q1(-1);
    q1=-q1;
    EXPECT_EQ(1,q1.U1());
    EXPECT_EQ(PRT_EVEN,q1.prt()); // No effects on the parity/fermionic parity.
    
    Qnum q2(PRTF_EVEN, 3, PRT_ODD);
    q2=-q2;
    EXPECT_EQ(-3,q2.U1());
    EXPECT_EQ(PRT_ODD,q2.prt()); //No effects on the parity/fermionic parity.
    EXPECT_EQ(PRTF_EVEN,q2.prtF()); //No effects on the parity/fermionic parity.
}

TEST(Qnum, Operation_LessThan){
    Qnum q1(3), q2(-1);
    
    EXPECT_TRUE(q2<q1);
    q1.assign(PRTF_ODD);
    q2.assign(PRTF_EVEN);
    
    EXPECT_TRUE(q2<q1);
    q1.assign(PRT_ODD);
    q2.assign(PRT_EVEN);
    
    EXPECT_TRUE(q2<q1);
}

TEST(Qnum, Operation_LessEq){
    Qnum q1(3), q2(-1);
    
    EXPECT_TRUE(q2<=q1);
    q1.assign(-1);
    EXPECT_TRUE(q2<=q1);
}

TEST(Qnum, Operation_Eq){
    Qnum q1(3), q2(2);
    
    EXPECT_FALSE(q2==q1);
    
    q1.assign(0,PRT_EVEN);
    q2.assign(0,PRT_ODD);
    EXPECT_FALSE(q2==q1);
    
    q1.assign(PRTF_EVEN,0,PRT_EVEN);
    q2.assign(PRTF_ODD,0,PRT_EVEN);
    EXPECT_FALSE(q2==q1);
}


TEST(Qnum, Operation_Mul){
    Qnum q1(PRTF_EVEN,3,PRT_ODD), q2(PRTF_ODD,2,PRT_ODD),q3;
    q3=q1*q2;
    EXPECT_EQ(5,q3.U1());
    EXPECT_EQ(PRT_EVEN,q3.prt());
    EXPECT_EQ(PRTF_ODD,q3.prtF());
    Qnum q4(-1,PRT_ODD), q5(2,PRT_ODD);
    
    q3=q4*q5;
    EXPECT_EQ(1,q3.U1());
    EXPECT_EQ(PRT_EVEN,q3.prt());
    EXPECT_NE(PRTF_ODD,q3.prtF());
    
}

