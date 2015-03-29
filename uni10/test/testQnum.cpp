//
//  Test_Qnum.cpp
//  
//
//  Created by Ying-Jer Kao on 3/20/15.
//
//

#include <gtest/gtest.h>
#include <iostream>
//#include <assert.h>
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

TEST(Qnum,CopyConstructor){
    // Qnum(const Qnum& _q);
    Qnum q1(1);
    
    Qnum q2(q1);
    
    ASSERT_EQ(q1,q2);
    
}


//int main(int argc, char **argv) {
//    InitGoogleTest(&argc, argv);
//    
//    return RUN_ALL_TESTS();
//}

