//
//  Test_Qnum.cpp
//
//
//  Created by Pochung Chen on 4/2/15.
//
//

#include <gtest/gtest.h>
#include <iostream>
//#include <EXPECT.h>
#include <map>
#include "uni10.hpp"
#include <time.h>
#include <vector>
using namespace uni10;
TEST(Bond,DefaultConstructor){
    // brief Default constuctor
    // Bond(bondType _type, size_t dim) : m_type(_type){

    Bond bd;

}

TEST(Bond,ConstructorTypeDim){
    // Create a Bond of type \c tp and dimension \c dim
    Bond bd(BD_IN, 100);

    EXPECT_EQ(BD_IN, bd.type());
    EXPECT_EQ(100,bd.dim());
}

TEST(Bond, ConstructorTypeQnums){
    // Bond(bondType tp, const std::vector<Qnum>& qnums);

    Qnum q1(1);
	Qnum q0(0);
	Qnum q_1(-1);
	// Create an array of Qnums for the states of a bond.
	std::vector<uni10::Qnum> qnums;
	qnums.push_back(q1);
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);

    // Constrcut Bond with Qnum array
	Bond bd(uni10::BD_OUT, qnums);

    // test dim
    EXPECT_EQ(6,bd.dim());

    // test Qlist
    std::vector<uni10::Qnum> qlist = bd.Qlist();
    EXPECT_EQ(q1,qlist[0]);
    EXPECT_EQ(q1,qlist[1]);
    EXPECT_EQ(q0,qlist[2]);
    EXPECT_EQ(q0,qlist[3]);
    EXPECT_EQ(q0,qlist[4]);
    EXPECT_EQ(q_1,qlist[5]);

    // test degeneracy
    std::map<Qnum, int> degs = bd.degeneracy();

    std::map<Qnum,int>::const_iterator it=degs.begin();
    EXPECT_EQ(q_1, it->first);
    EXPECT_EQ(1, it->second);

    ++it;
    EXPECT_EQ(q0,it->first);
    EXPECT_EQ(3,it->second);

    ++it;
    EXPECT_EQ(q1,it->first);
    EXPECT_EQ(2,it->second);

}

TEST(Bond, CopyConstructor){
    // Bond(const Bond& bd);
    Bond bd1(BD_IN, 100);
    Bond bd2(bd1);

    EXPECT_EQ(bd1,bd2);
}

TEST(Bond, ChangeBondType){
    // Bond& change(bondType tp);

    Bond bd(BD_IN, 100);
    bd.change(BD_OUT);
    EXPECT_EQ(BD_OUT,bd.type());
    // test if the qnum is inverted
}

TEST(Bond, DummyChangeBondType){
    // Bond& dummy_change(bondType tp);
}

TEST(Bond, combine){
    // Bond& combine(Bond bd);
    uni10::Qnum q2(2);
    uni10::Qnum q1(1);
	uni10::Qnum q0(0);
    uni10::Qnum q_1(-1);
	uni10::Qnum q_2(-2);
	// Create an array of Qnums for the states of a bond.
	std::vector<uni10::Qnum> qnums;
	qnums.push_back(q1);
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);

	// Constrcut first bond
	uni10::Bond bd(uni10::BD_IN, qnums);

    // Construct another bond
	qnums.clear();
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);
	uni10::Bond bd2(uni10::BD_IN, qnums);
    bd2.combine(bd);

  //  std::cout<<"Degeneracies of bd2 after combining bd: "<<std::endl;
    std::map<uni10::Qnum, int> degs;
	degs = bd2.degeneracy();
	//for(std::map<uni10::Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it)
	//	std::cout<<it->first<<": "<<it->second<<std::endl;
	//std::cout<<std::endl;

    // test bond type
    EXPECT_EQ(BD_IN,bd2.type());
    // test bond dimension
    EXPECT_EQ(24, bd2.dim());
    // test degeneracy
    std::map<Qnum,int>::const_iterator it=degs.begin();
    EXPECT_EQ(q_2, it->first);
    EXPECT_EQ(1, it->second);

    ++it;
    EXPECT_EQ(q_1,it->first);
    EXPECT_EQ(5,it->second);

    ++it;
    EXPECT_EQ(q0,it->first);
    EXPECT_EQ(9,it->second);

    ++it;
    EXPECT_EQ(q1,it->first);
    EXPECT_EQ(7,it->second);

    ++it;
    EXPECT_EQ(q2,it->first);
    EXPECT_EQ(2,it->second);
}
