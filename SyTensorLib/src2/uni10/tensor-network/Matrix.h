#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <string.h>
#include <math.h>
//Type of Matrix
//#include <uni10/tensor-network/UniTensor.h>
namespace uni10{

class Matrix_t {
	public:
		Matrix_t(int _Rnum, int _Cnum, bool _diag=false);
		Matrix_t(int _Rnum, int _Cnum, double* _elem, bool _diag=false);
		Matrix_t(const Matrix_t& _m);
		~Matrix_t();
		int row()const;
		int col()const;
		bool isDiag()const{return diag;};
		size_t getElemNum()const;
		Matrix_t& operator=(const Matrix_t& _m);
		friend Matrix_t operator* (const Matrix_t& Ma, const Matrix_t& Mb);
		Matrix_t& operator*= (const Matrix_t& Mb);
		friend std::ostream& operator<< (std::ostream& os, const Matrix_t& b);
		std::vector<Matrix_t> diagonalize();
		std::vector<Matrix_t> svd();
		void orthoRand();
		void bzero();
		void transpose();
		double trace();
		void save(const std::string fname);
		void load(const std::string fname);
		friend Matrix_t operator*(const Matrix_t& Ma, double a);
		friend Matrix_t operator*(double a, const Matrix_t& Ma){return Ma * a;};
		friend bool operator== (const Matrix_t& m1, const Matrix_t& m2);
		void operator*= (double a);
		friend Matrix_t operator+(const Matrix_t& Ma, const Matrix_t& Mb);
		void operator+= (const Matrix_t& Mb);
		double& operator[](size_t idx);
		double* elem;
	private:
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		size_t elemNum;
		bool diag;
};

};	/* namespace uni10 */	
#endif /* MATRIX_H */
