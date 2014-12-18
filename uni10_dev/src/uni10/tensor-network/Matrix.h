/****************************************************************************
*  @file Matrix.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    Yun-Da Hsieh, Pochung Chen and Ying-Jer Kao
*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Header file for Matrix class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <string.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <uni10/data-structure/Block.h>
//Type of Matrix
namespace uni10{

class Matrix: public Block {
	public:
		Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=true);
		Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag=false, bool src_ongpu=false);
		Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag=false, bool src_ongpu=false);
		Matrix(const Matrix& _m);
		Matrix();
		~Matrix();
		Matrix& operator=(const Matrix& _m);
		//size_t row()const;
		//size_t col()const;
		//bool isDiag()const{return diag;};
		//bool isOngpu()const{return ongpu;};
		//size_t elemNum()const;
		double& operator[](size_t idx);
		double& at(size_t i, size_t j);
		//double* getElem()const;
		double* getHostElem();
		void setElem(const double* elem, bool _ongpu = false);
		void setElem(const std::vector<double>& elem, bool _ongpu = false);
		Matrix& resize(size_t row, size_t col);
		//void save(const std::string& fname)const;
		void load(const std::string& fname);
		void identity();
		void set_zero();
		void randomize();
		void orthoRand();
		Matrix& transpose();
		//std::vector<Matrix> eigh()const;
		std::vector<Matrix> svd()const;
    size_t lanczosEigh(double& E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15)const;
		double trace()const;
		double norm()const;
		double sum()const;
		Matrix& operator*= (double a);
		Matrix& operator*= (const Matrix& Mb);
		Matrix& operator+= (const Matrix& Mb);
		friend Matrix takeExp(double a, const Matrix& mat);
		friend Matrix operator* (const Matrix& Ma, const Matrix& Mb);
		friend Matrix operator*(const Matrix& Ma, double a);
		friend Matrix operator*(double a, const Matrix& Ma){return Ma * a;};
		friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);
		friend bool operator== (const Matrix& m1, const Matrix& m2);
		friend std::ostream& operator<< (std::ostream& os, const Matrix& b);
		bool toGPU();
	private:
		void init(bool togpu);
		void init(const double* elem, bool _ongpu);
    /*
		size_t Rnum;		//number of rows of the block
		size_t Cnum;		//number of columns of the block
		double* m_elem;
		size_t m_elemNum;
		bool diag;
		bool ongpu;
    */
};
Matrix takeExp(double a, const Matrix& mat);
Matrix otimes(const Matrix& Ma, const Matrix& Mb);
};	/* namespace uni10 */
#endif /* MATRIX_H */
