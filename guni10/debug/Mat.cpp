#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
	Qnum q0;
	Bond bdi(BD_IN, 5);
	Bond bdo(BD_OUT, 5);
	vector<Bond> bonds(2, bdi);
	bonds.push_back(bdo);
	bonds.push_back(bdo);

	UniTensor T(bonds);
	T.randomize();
	Matrix D(12, 20, false, false);
	D.load("gpu_elem");
	vector<Matrix>rets = D.svd();
	cout<<D;
	
	Matrix S(12, 12, false, true);
	S.load("square_elem");
	cout<<S;

	cout<<"TRACE = "<<S.trace()<<endl;
	
	Matrix D1(12, 12, false, true);	
	Matrix D2(12, 12, false, false);	
	//DD.randomize();
	D1.load("square_elem");
	D2.load("square_elem");
	//cout<<D1;
	cout<<"TRACE = "<<D1.trace()<<endl;
	Matrix D3 = D1 + D2;	
	//cout<<D3;

	//D3.getHostElem();
	//cout<<"D3.ongpu = "<<D3.isOngpu()<<", D3[0] = "<<D3[0]<<endl;

	
	cout<<"ongpu = "<<D1.isOngpu()<<endl;
	rets = D1.diagonalize();
	cout<<"------------------------\n";
	cout<<rets[0];
	cout<<rets[1];
}
