#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q_10(-1, 0);
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
						   0,      0,      0, 1.0/4
					  };
	SyTensor_t H0(bonds, "H0");
	H0.addRawElem(H_elem);
	vector<Qnum_t> H0_qnums = H0.qnums();
	for(int q = 0; q < H0_qnums.size(); q++){
		Matrix_t mat = H0.getBlock(H0_qnums[q]);
		cout<<mat;
	}
	Qnum_t q00(0, 0);
	Matrix_t H00 = H0.getBlock(q00);
	cout <<H00;
	H00.elem[2] = 9;
	H0.putBlock(q00, H00);
	cout<< H0;



	//----------
	Matrix_t mat(4, 4, H_elem);	
	cout << mat;
	Matrix_t mat1 = mat;

	mat1 *= mat;
	cout<< mat1;

	vector<Matrix_t> outs = mat.diagonalize();
	cout << mat;

	//-------
		
	Matrix_t Msvd(4, 5);
	Msvd.orthoRand();
	cout << Msvd;
	cout << "===== Transpose =====\n";
	Msvd.transpose();
	cout << Msvd;
	cout << "=====================\n";
	vector<Matrix_t> outs2 = Msvd.svd();
	cout << outs2[0];
	cout << outs2[1];
	cout << outs2[2];
	Matrix_t ret = outs2[0] * outs2[2];
	//ret *= outs2[2];
	cout<<ret;
}
