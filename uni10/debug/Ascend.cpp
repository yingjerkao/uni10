#include <iostream>
#include <assert.h>
#include <vector>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
//#define FERMION 1


int main(){
	Qnum q10(1, PRT_EVEN);
	Qnum q_10(-1, PRT_EVEN);
	Qnum q30(3, PRT_EVEN);
#ifdef FERMION
	Qnum q_11(PRTF_ODD, -1, PRT_EVEN);
	Qnum q11(PRTF_ODD, 1, PRT_EVEN);
	Qnum q_31(PRTF_ODD, -3, PRT_EVEN);
#else
	Qnum q_11(-1, PRT_ODD);
	Qnum q11(1, PRT_ODD);
	Qnum q_31(-3, PRT_ODD);
#endif
	vector<Bond> bonds;
	vector<Qnum> qnums;
	qnums.push_back(q10);qnums.push_back(q_11);
	Bond bdr(BD_IN, qnums);
	Bond bdc(BD_OUT, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	UniTensor H0(bonds, "Ob");
	H0.setRawElem(H_elem);

	UniTensor U(bonds, "U");
	U.orthoRand();

	bonds.clear();
	vector<Qnum> qnums1;
	qnums1.push_back(q30); qnums1.push_back(q11); qnums1.push_back(q11); qnums1.push_back(q11);
	qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_10); qnums1.push_back(q_31);
	Bond bdr1(BD_IN, qnums1);
	Bond bdc1(BD_OUT, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	UniTensor W1(bonds, "W1");
	W1.orthoRand();
	UniTensor W2(bonds, "W2");
	W2.orthoRand();

	UniTensor H1;
	CNetwork asdL("AscendL");
  asdL.putTensor("W1", W1);
  asdL.putTensorT("W1T", W1);
  asdL.putTensor("W2", W2);
  asdL.putTensorT("W2T", W2);
  asdL.putTensor(2, U);
  asdL.putTensorT("UT", U);
  asdL.putTensor("Ob", H0);
  cout<<asdL.profile();
	H1 = asdL.launch();
	H0.profile();
  cout<<H1;
  Qnum q01(PRTF_ODD, 0);
  /*
  Matrix blk = H1.getBlock(q01);
  for(map<Qnum, Block>::const_iterator it = H1.const_getBlocks().begin(); it != H1.const_getBlocks().end(); it++)
    cout<<it->second.trace()<<endl;
  cout<<H1.trace()<<endl;
  Matrix mat = blk;
  */

}

