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
    ASSERT_EQ(q.U1(),0);
    ASSERT_EQ(q.prt(),PRT_EVEN);
    ASSERT_FALSE(q.isFermionic());
}



TEST(Qnum,DefaultFermionConstructor){
    // Qnum(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);
    Qnum q(PRTF_ODD);
    ASSERT_EQ(q.prtF(),PRTF_ODD);
    ASSERT_EQ(q.U1(),0);
    ASSERT_EQ(q.prt(),PRT_EVEN);
    ASSERT_TRUE(q.isFermionic());
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
    
    ASSERT_EQ(q1,q2);
    
}
TEST(Qnum,isFermionic){
    // Qnum(int _U1=0, parityType _prt=PRT_EVEN);
    Qnum q;
  
    ASSERT_EQ(q.U1(),0);
    ASSERT_EQ(q.prt(),PRT_EVEN);
// PRTF_ODD is defined, so fermionic is set true
    ASSERT_TRUE(q.isFermionic()); 
}

TEST(Qnum,assign){
    Qnum q;
    // assign(int _U1 = 0, parityType _prt = PRT_EVEN);
    q.assign(3);
    ASSERT_EQ(q.U1(),3);
    ASSERT_EQ(q.prt(),PRT_EVEN);
    
    q.assign(0,PRT_ODD);
    ASSERT_EQ(q.prt(),PRT_ODD);
    ASSERT_EQ(q.U1(),0);
    
    q.assign(-1);
    ASSERT_EQ(q.prt(),PRT_EVEN); // Default value PRT_EVEN assigned
    ASSERT_EQ(q.U1(),-1);
    
    Qnum qf;
    //assign(parityFType _prtF, int _U1 = 0, parityType _prt = PRT_EVEN);
    
    qf.assign(PRTF_EVEN, 1, PRT_ODD);
    
    ASSERT_EQ(qf.prtF(),PRTF_EVEN);
    ASSERT_EQ(qf.U1(),1);
    ASSERT_EQ(qf.prt(), PRT_ODD);

    qf.assign(PRTF_ODD);
    
    ASSERT_EQ(qf.prtF(), PRTF_ODD);
    ASSERT_EQ(qf.prt(),PRT_EVEN); // Default values are also implicitly assigned
    ASSERT_EQ(qf.U1(),0);

}

TEST(Qnum,hash){
    Qnum q(PRTF_ODD,100,PRT_EVEN);
    ASSERT_EQ(q.hash(), 401);
}


TEST(Qnum, Operation_Negation){
    Qnum q1(-1);
    q1=-q1;
    ASSERT_EQ(q1.U1(),1);
    ASSERT_EQ(q1.prt(),PRT_EVEN); // No effects on the parity/fermionic parity.
    
    Qnum q2(PRTF_EVEN, 3, PRT_ODD);
    q2=-q2;
    ASSERT_EQ(q2.U1(),-3);
    ASSERT_EQ(q2.prt(),PRT_ODD); //No effects on the parity/fermionic parity.
    ASSERT_EQ(q2.prtF(),PRTF_EVEN); //No effects on the parity/fermionic parity.
}

TEST(Qnum, Operation_LessThan){
    Qnum q1(3), q2(-1);
    
    ASSERT_TRUE(q2<q1);
    q1.assign(PRTF_ODD);
    q2.assign(PRTF_EVEN);
    
    ASSERT_TRUE(q2<q1);
    q1.assign(PRT_ODD);
    q2.assign(PRT_EVEN);
    
    ASSERT_TRUE(q2<q1);
}

TEST(Qnum, Operation_LessEq){
    Qnum q1(3), q2(-1);
    
    ASSERT_TRUE(q2<=q1);
    q1.assign(-1);
    ASSERT_TRUE(q2<=q1);
}

TEST(Qnum, Operation_Eq){
    Qnum q1(3), q2(2);
    
    ASSERT_FALSE(q2==q1);
    
    q1.assign(0,PRT_EVEN);
    q2.assign(0,PRT_ODD);
    ASSERT_FALSE(q2==q1);
    
    q1.assign(PRTF_EVEN,0,PRT_EVEN);
    q2.assign(PRTF_ODD,0,PRT_EVEN);
    ASSERT_FALSE(q2==q1);
}


TEST(Qnum, Operation_Mul){
    Qnum q1(PRTF_EVEN,3,PRT_ODD), q2(PRTF_ODD,2,PRT_ODD),q3;
    q3=q1*q2;
    ASSERT_EQ(q3.U1(),5);
    ASSERT_EQ(q3.prt(),PRT_EVEN);
    ASSERT_EQ(q3.prtF(),PRTF_ODD);
    Qnum q4(-1,PRT_ODD), q5(2,PRT_ODD);
    
    q3=q4*q5;
    ASSERT_EQ(q3.U1(),1);
    ASSERT_EQ(q3.prt(),PRT_EVEN);
    ASSERT_NE(q3.prtF(),PRTF_ODD);
    
}

