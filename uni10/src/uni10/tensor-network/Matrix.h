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

class Matrix {
	public:
		Matrix(int _Rnum, int _Cnum, bool _diag=false);
		Matrix(int _Rnum, int _Cnum, double* _elem, bool _diag=false);
		Matrix(const Matrix& _m);
		~Matrix();
		int row()const;
		int col()const;
		bool isDiag()const{return diag;};
		size_t elemNum()const;
		Matrix& operator=(const Matrix& _m);
		friend Matrix operator* (const Matrix& Ma, const Matrix& Mb);
		Matrix& operator*= (const Matrix& Mb);
		friend std::ostream& operator<< (std::ostream& os, const Matrix& b);
		std::vector<Matrix> diagonalize();
		std::vector<Matrix> svd();
		void addElem(double* elem);
		void randomize();
		void orthoRand();
		void set_zero();
		void transpose();
		double trace();
		void save(const std::string& fname);
		void load(const std::string& fname);
		friend Matrix operator*(const Matrix& Ma, double a);
		friend Matrix operator*(double a, const Matrix& Ma){return Ma * a;};
		friend bool operator== (const Matrix& m1, const Matrix& m2);
		Matrix& operator*= (double a);
		friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);
		Matrix& operator+= (const Matrix& Mb);
		double& operator[](size_t idx);
		double* elem()const;
		double& at(int i, int j);
	private:
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		double* m_elem;
		size_t m_elemNum;
		bool diag;
};

};	/* namespace uni10 */	
#endif /* MATRIX_H */
