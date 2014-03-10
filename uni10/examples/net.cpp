#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q_10(-1, 0);
	Qnum_t q30(3, 0);
	Qnum_t q_30(-3, 0);
	vector<SyTensor_t*> SyTptrs;
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q10);qnums.push_back(q_10);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};
	int label_H0[] = {3, 7, 4, 8};
	SyTensor_t H0(bonds, label_H0, "H0");
	H0.addRawElem(H_elem);
	cout << H0;

	SyTensor_t U(bonds, "U");
	U.orthoRand();
	SyTensor_t UT = U;
	UT.transpose();
	//int label_U[] = {2, 6, 3, 7};
	//U.addLabel(label_U);
	//int label_UT[] = {4, 8, 5, 9};
	//UT.addLabel(label_UT);
	//UT.setName("UT");
	
	bonds.clear();
	vector<Qnum_t> qnums1;
	qnums1.push_back(q30);qnums1.push_back(q10);qnums1.push_back(q10);qnums1.push_back(q10);
	qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_10);qnums1.push_back(q_30);
	Bond_t bdr1(BD_ROW, qnums1);
	bonds.clear();
	bonds.push_back(bdr1);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	SyTensor_t W1(bonds, "W1");
	W1.orthoRand();
	cout << W1;
	SyTensor_t W1T = W1;
	W1T.transpose();
	int label_W1[] = {-1, 0, 1, 2};
	W1.addLabel(label_W1);
	int label_W1T[] = {0, 1, 5, -3};
	W1T.addLabel(label_W1T);
	//W1T.setName("W1T");

	SyTensor_t W2(bonds);
	W2.orthoRand();
	SyTensor_t W2T = W2;
	W2T.transpose();
	int label_W2[] = {-2, 6, 10, 11};
	W2.addLabel(label_W2);
	int label_W2T[] = {9, 10, 11, -4};
	W2T.addLabel(label_W2T);
	//W2T.setName("W2T");
	int label_out[] = {-1, -2, -3, -4};

	/*
	SyTensor_t H1 = W1 * W1T;
	SyTensor_t tmp = UT * U;
	tmp *= W2T;
	tmp *= H0;
	tmp *= W2;
	H1 *= tmp;	
	H1.reshape(label_out, 2);
	cout<<H1;
	*/


	/*
	printRawElem(W1);
	printRawElem(W2);
	Network_t net;
	net.add(U);
	net.add(H0);
	net.add(UT);
	net.add(W1);
	net.add(W1T);
	net.add(W2);
	net.add(W2T);
	SyTensor_t Tret = net.launch(label_out, 2);
	Tret.save("H1");
	*/

	Network_t net1("AscendC");
	net1.replaceWith(0, &W1);
	net1.replaceWith(1, &W1T);
	net1.replaceWith(2, &U);
	net1.replaceWith(3, &H0);
	net1.replaceWith(4, &UT);
	net1.replaceWith(5, &W2);
	net1.replaceWith(6, &W2T);
	SyTensor_t Tret1 = net1.launch();
	Tret1.save("H11");
	cout<<net1;

	//cout<<Tret;
	//cout<<Tret;
	//printRawElem(Tret);
	W1.check();
	cout<<Tret1;}
