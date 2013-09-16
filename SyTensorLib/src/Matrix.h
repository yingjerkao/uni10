#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <string.h>
#include <map>
using namespace std;
//Type of Matrix
#include "myLapack.h"

class SyTensor_t;
class Matrix_t {
	public:
		Matrix_t(int _Rnum, int _Cnum, bool _diag=false);
		Matrix_t(int _Rnum, int _Cnum, double* _elem, bool _diag=false);
		Matrix_t(const Matrix_t& _m);
		~Matrix_t();
		int row();
		int col();
		bool isDiag(){return diag;};
		Matrix_t& operator=(const Matrix_t& _m);
		friend Matrix_t operator* (const Matrix_t& Ma, const Matrix_t& Mb);
		void operator*= (const Matrix_t& Mb);
		friend ostream& operator<< (ostream& os, const Matrix_t& b);
		vector<Matrix_t> diagonalize();
		vector<Matrix_t> svd();
		void orthoRand();
		void transpose();
		friend Matrix_t operator*(const Matrix_t& Ma, double a);
		friend Matrix_t operator*(double a, const Matrix_t& Ma){return Ma * a;};
		void operator*= (double a);
		friend Matrix_t operator+(const Matrix_t& Ma, const Matrix_t& Mb);
		void operator+= (const Matrix_t& Mb);
		double* elem;
	private:
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		int elemNum;
		bool diag;
};
