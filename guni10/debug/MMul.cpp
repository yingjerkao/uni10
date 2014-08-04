#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
	Matrix M1(12, 12, false, false);	
	Matrix M2(12, 16, false, false);	
	M1.randomize();
	M2.randomize();
	cout<<M1 * M2;
	Matrix Ma(M1);
	Matrix Mb(M2);
	Matrix Mc = Ma * Mb;
	cout<<Mc;
	cout<<"isOngpu = "<<Mc.isOngpu()<<endl;

	double H_elem[] = {\
		1.0/4,      0,      0,     0,
		    0, -1.0/4,  1.0/2,     0,
	        0,  1.0/2, -1.0/4,     0,
	        0,      0,      0, 1.0/4};
	Matrix H(4, 4, H_elem);
	cout<<H;
	vector<Matrix> rets = H.diagonalize();
	cout<<rets[0];
	cout<<rets[1];
	H.getHostElem();
	cout<<H;
	cout<<takeExp(-0.01, H);

	Matrix MS(12, 16, false, true);
	MS.randomize();
	cout<<MS;
	MS.resize(8, 16);
	cout<<MS;
	MS.resize(10, 8);
	cout<<MS;
}
