/****************************************************************************
*  @file CMakeLists.txt
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
*  @brief Main specification file for CMake
*  @author Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <uni10/tensor-network/Matrix.h>
#include <uni10/numeric/uni10_lapack.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10{
std::ostream& operator<< (std::ostream& os, const Matrix& m){
	os << m.Rnum << " x " << m.Cnum << " = " << m.m_elemNum;
	if(m.diag)
		os << ", Diagonal";
	os <<std::endl << std::endl;
	for(int i = 0; i < m.Rnum; i++){
		for(int j = 0; j < m.Cnum; j++)
			if(m.diag){
				if(i == j)
					os << std::setw(7) << std::fixed << std::setprecision(3) << m.m_elem[i];
				else
					os << std::setw(7) << std::fixed << std::setprecision(3) << 0.0;
			}
			else
				os << std::setw(7) << std::fixed << std::setprecision(3) << m.m_elem[i * m.Cnum + j];
		os << std::endl << std::endl;
	}
	return os;
}

Matrix::Matrix(): Rnum(0), Cnum(0), m_elemNum(0), diag(false), m_elem(NULL){
}
Matrix::Matrix(const Matrix& _m): Rnum(_m.Rnum), Cnum(_m.Cnum), m_elemNum(_m.m_elemNum), diag(_m.diag), m_elem(NULL){
	if(m_elemNum){
		m_elem = (double*)malloc(m_elemNum * sizeof(double));
		memcpy(m_elem, _m.m_elem, m_elemNum * sizeof(double));
	}
}

Matrix::Matrix(int _Rnum, int _Cnum, double* _elem, bool _diag): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), m_elem(NULL){
	if(_diag)
		m_elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	if(m_elemNum){
		m_elem = (double*)malloc(m_elemNum * sizeof(double));
		memcpy(m_elem, _elem, m_elemNum * sizeof(double));
	}
}

Matrix::Matrix(int _Rnum, int _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), m_elemNum(_Rnum * _Cnum), diag(_diag), m_elem(NULL){
	if(_diag)
		m_elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	if(m_elemNum){
		m_elem = (double*)malloc(m_elemNum * sizeof(double));
		memset(m_elem, 0, m_elemNum * sizeof(double));
	}
}

Matrix& Matrix::operator=(const Matrix& _m){
	Rnum = _m.Rnum;
	Cnum = _m.Cnum;
	m_elemNum = _m.m_elemNum;
	diag = _m.diag;
	m_elem = (double*)realloc(m_elem, Rnum * Cnum * sizeof(double));
	memcpy(m_elem, _m.m_elem, Rnum * Cnum * sizeof(double));
	return *this;
}

Matrix::~Matrix(){
	if(m_elem != NULL)
		free(m_elem);
}

int Matrix::row()const{
	return Rnum;
}

int Matrix::col()const{
	return Cnum;
}
size_t Matrix::elemNum()const{
	return m_elemNum;
}

Matrix operator* (const Matrix& Ma, const Matrix& Mb){
	assert(Ma.Cnum == Mb.Rnum);
	if((!Ma.diag) && (!Mb.diag)){
		Matrix Mc(Ma.Rnum, Mb.Cnum);
		myDgemm(Ma.m_elem, Mb.m_elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.m_elem);
		return Mc;
	}
	else if(Ma.diag && (!Mb.diag)){
		Matrix Mc(Ma.Rnum, Mb.Cnum);
		for(int i = 0; i < Ma.m_elemNum; i++)
			for(int j = 0; j < Mb.Cnum; j++)
				Mc.m_elem[i * Mb.Cnum + j] = Ma.m_elem[i] * Mb.m_elem[i * Mb.Cnum + j];
		return Mc;
	}
	else if((!Ma.diag) && Mb.diag){
		Matrix Mc(Ma.Rnum, Mb.Cnum);
		for(int i = 0; i < Ma.Rnum; i++)
			for(int j = 0; j < Mb.m_elemNum; j++)
				Mc.m_elem[i * Mb.Cnum + j] = Ma.m_elem[i * Ma.Cnum + j] * Mb.m_elem[j];
		return Mc;
	}
	else{
		Matrix Mc(Ma.Rnum, Mb.Cnum, true);
		for(int i = 0; i < Ma.Rnum; i++)
			Mc.m_elem[i] = Ma.m_elem[i] * Mb.m_elem[i];
		return Mc;
	}
}
bool operator== (const Matrix& m1, const Matrix& m2){
	double diff;
	if(m1.m_elemNum == m2.m_elemNum){
		for(int i = 0; i < m1.m_elemNum; i++){
			diff = fabs(m1.m_elem[i] - m2.m_elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}


Matrix& Matrix::operator*= (const Matrix& Mb){
	return *this = *this * Mb;
}

void Matrix::setElem(double* elem){
	memcpy(m_elem, elem, m_elemNum * sizeof(double));
}

std::vector<Matrix> Matrix::diagonalize()const{
	assert(Rnum == Cnum);
	assert(!diag);
	std::vector<Matrix> outs;
	Matrix Eig(Rnum, Cnum, true);
	Matrix EigV(Rnum, Cnum);
	syDiag(m_elem, Rnum, Eig.m_elem, EigV.m_elem);
	outs.push_back(Eig);
	outs.push_back(EigV);
	return outs;
}

std::vector<Matrix> Matrix::svd() const{
	assert(!diag);
	std::vector<Matrix> outs;
	int min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
	Matrix U(Rnum, min);
	Matrix S(min, min, true);
	Matrix VT(min, Cnum);
	myDgesvd(m_elem, Rnum, Cnum, U.m_elem, S.m_elem, VT.m_elem);
	outs.push_back(U);
	outs.push_back(S);
	outs.push_back(VT);
	return outs;
}

void Matrix::randomize(){
	randomNums(m_elem, m_elemNum, 0);
}
void Matrix::orthoRand(){
	if(!diag){
		if(Rnum <= Cnum)
			orthoRandomize(m_elem, Rnum, Cnum);
		else{
			Matrix M(Cnum, Rnum);
			orthoRandomize(M.getElem(), Cnum, Rnum);
			myTranspose(M.getElem(), Cnum, Rnum, m_elem, 0);
		}

	}
}

void Matrix::set_zero(){
	if(m_elemNum)
		memset(m_elem, 0, m_elemNum * sizeof(double));
}

Matrix operator*(const Matrix& Ma, double a){
	Matrix Mb(Ma);
	vecScal(a, Mb.m_elem, Mb.m_elemNum);
	return Mb;
}

Matrix& Matrix::operator*= (double a){
	vecScal(a, m_elem, m_elemNum);
	return *this;
}

Matrix operator+(const Matrix& Ma, const Matrix& Mb){
	Matrix Mc(Ma);
	vecAdd(Mb.m_elem, Mc.m_elem, Ma.m_elemNum);
	return Mc;
}

Matrix& Matrix::operator+= (const Matrix& Mb){
	vecAdd(Mb.m_elem, m_elem, m_elemNum);
	return *this;
}

void Matrix::transpose(){
	if(!diag){
		double* oldElem = (double*)malloc(m_elemNum * sizeof(double));
		memcpy(oldElem, m_elem, m_elemNum * sizeof(double));
		myTranspose(oldElem, Rnum, Cnum, m_elem, 0);
		free(oldElem);
	}
	int tmp = Rnum;
	Rnum = Cnum;
	Cnum = tmp;
}
double Matrix::norm(){
	double nm = 0;
	for(int i = 0; i < m_elemNum; i++)
		nm += m_elem[i] * m_elem[i];
	return sqrt(nm);
}
double Matrix::sum(){
	double sm = 0;
	for(int i = 0; i < m_elemNum; i++)
		sm += m_elem[i];
	return sm;
}
double Matrix::trace(){
	assert(Rnum == Cnum);
	double sum = 0;
	if(diag)
		for(int i = 0; i < m_elemNum; i++)
			sum += m_elem[i];
	else
		for(int i = 0; i < Rnum; i++)
			sum += m_elem[i * Cnum + i];
	return sum;
}
void Matrix::save(const std::string& fname){
	FILE *fp = fopen(fname.c_str(), "w");
	assert(fp != NULL);
	fwrite(m_elem, sizeof(double), m_elemNum, fp);
	fclose(fp);
}
void Matrix::load(const std::string& fname){
	FILE *fp = fopen(fname.c_str(), "r");
	assert(fp != NULL);
	fread(m_elem, sizeof(double), m_elemNum, fp);
	fclose(fp);
}
double& Matrix::operator[](size_t idx){
	assert(idx < m_elemNum);
	return m_elem[idx];
}
double* Matrix::getElem()const{
	return m_elem;
}
double& Matrix::at(int r, int c){
	assert(r < Rnum);
	assert(c < Cnum);
	if(diag){
		assert(r == c && r < m_elemNum);
		return m_elem[r];
	}
	else
		return m_elem[r * Cnum + c];
}

Matrix takeExp(double a, const Matrix& mat){
	std::vector<Matrix> rets = mat.diagonalize();
	Matrix UT = rets[1];
	UT.transpose();
	for(int i = 0; i < rets[0].row(); i++)
		rets[0][i] = exp(a * rets[0][i]);
	return UT * rets[0] * rets[1];
}
};	/* namespace uni10 */
