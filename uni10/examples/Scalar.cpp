#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
//using namespace uni10;

int main(){
	uni10::Qnum q1(1);
	uni10::Qnum q0(0);
	uni10::Qnum q_1(-1);
	uni10::UniTensor T(9);
	T.printRawElem();
	vector<uni10::Qnum> qnums;
	//qnums.push_back(q1);
	qnums.push_back(q1);
	qnums.push_back(q_1);
	//qnums.push_back(q_1);
	uni10::Bond BI(uni10::BD_IN, qnums);
	uni10::Bond BO(uni10::BD_OUT, qnums);
	vector<uni10::Bond> bonds;
	bonds.push_back(BI);
	bonds.push_back(BI);
	bonds.push_back(BO);
	bonds.push_back(BO);

	uni10::UniTensor T1(bonds);
	uni10::UniTensor T2(bonds);
	T1.randomize();
	T2.randomize();
	//T1.save("Data/Tata/T1.ut");
	uni10::UniTensor Tf1 = T1 * T2;
	int labels[] = {0, 1, 4, 5};
	int per_labels[] = {2, 3, 0, 1};
	T2.addLabel(labels);
	//T1.permute(per_labels, 2);
	T1.transpose();
	T1.addLabel(per_labels);
	cout<<T1;
	T1.printRawElem();
	cout<<T2;
	T2.printRawElem();
	uni10::UniTensor Tf2 = T1 * T2;
	cout<<Tf2;
	Tf2.printRawElem();
	cout<<Tf1;
	cout<<Tf2.trace()<<endl;
	uni10::UniTensor T3("Data/Tata/T1.ut");
	cout<<T3;
}
