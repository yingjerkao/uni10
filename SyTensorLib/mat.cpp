#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "Matrix.h"

int main(){
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4
					  };
	Matrix_t mat(4, 4, H_elem);	
	cout << mat;
	Matrix_t mat1 = mat;

	mat1 *= mat;
	cout<< mat1;

	vector<Matrix_t> outs = mat.diagonalize();
	cout << mat;
	double D_elem[] = {1, 2, 3, 4};
	Matrix_t D(5, 4, D_elem, true);	
		
	Matrix_t Msvd(4, 5);
	Msvd.orthoRand();
	cout << Msvd;
	vector<Matrix_t> outs2 = Msvd.svd();
	cout << outs2[0];
	cout << outs2[1];
	cout << outs2[2];
	Matrix_t ret = outs2[0];
	ret *= outs2[1];
	cout<<ret;
	ret *= outs2[2];
	cout<<ret;
}
