/****************************************************************************
*  @file Matrix.cpp
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
*  @brief Implementation file for Matrix class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include "Matrix.h"

ostream& operator<< (ostream& os, const Matrix_t& m){
	os << endl << m.Rnum << " x " << m.Cnum << " = " << m.elemNum;
	if(m.diag)
		os << ", Diagonal";
	os <<endl << endl;
	for(int i = 0; i < m.Rnum; i++){
		for(int j = 0; j < m.Cnum; j++)
			if(m.diag)
				if(i == j)
					os << setw(7) << fixed << setprecision(3) << m.elem[i];
				else
					os << setw(7) << fixed << setprecision(3) << 0.0;
			else
				os << setw(7) << fixed << setprecision(3) << m.elem[i * m.Cnum + j];
		os << endl << endl;
	}
	return os;
}

Matrix_t::Matrix_t(const Matrix_t& _m): Rnum(_m.Rnum), Cnum(_m.Cnum), elemNum(_m.elemNum), diag(_m.diag){
	elem = (double*)malloc(elemNum * sizeof(double));
	memcpy(elem, _m.elem, elemNum * sizeof(double));
}

Matrix_t::Matrix_t(int _Rnum, int _Cnum, double* _elem, bool _diag): Rnum(_Rnum), Cnum(_Cnum), elemNum(_Rnum * _Cnum), diag(_diag), elem(NULL){
	if(_diag)
		elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	elem = (double*)malloc(elemNum * sizeof(double));
	memcpy(elem, _elem, elemNum * sizeof(double));
}

Matrix_t::Matrix_t(int _Rnum, int _Cnum, bool _diag): Rnum(_Rnum), Cnum(_Cnum), elemNum(_Rnum * _Cnum), diag(_diag), elem(NULL){
	if(_diag)
		elemNum = _Rnum < _Cnum ? _Rnum : _Cnum;
	elem = (double*)malloc(elemNum * sizeof(double));
	memset(elem, 0, elemNum * sizeof(double));
}

Matrix_t& Matrix_t::operator=(const Matrix_t& _m){
	Rnum = _m.Rnum;
	Cnum = _m.Cnum;
	elemNum = _m.elemNum;
	diag = _m.diag;
	elem = (double*)realloc(elem, Rnum * Cnum * sizeof(double));
	memcpy(elem, _m.elem, Rnum * Cnum * sizeof(double));
	return *this;
}

Matrix_t::~Matrix_t(){
	free(elem);
}

int Matrix_t::row()const{
	return Rnum;
}

int Matrix_t::col()const{
	return Cnum;
}

Matrix_t operator* (const Matrix_t& Ma, const Matrix_t& Mb){
	assert(Ma.Cnum == Mb.Rnum);
	if((!Ma.diag) && (!Mb.diag)){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		myDgemm(Ma.elem, Mb.elem, Ma.Rnum, Mb.Cnum, Ma.Cnum, Mc.elem);
		return Mc;
	}
	else if(Ma.diag && (!Mb.diag)){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		for(int i = 0; i < Ma.elemNum; i++)
			for(int j = 0; j < Mb.Cnum; j++)
				Mc.elem[i * Mb.Cnum + j] = Ma.elem[i] * Mb.elem[i * Mb.Cnum + j];
		return Mc;
	}
	else if((!Ma.diag) && Mb.diag){
		Matrix_t Mc(Ma.Rnum, Mb.Cnum);	
		for(int i = 0; i < Ma.Rnum; i++)
			for(int j = 0; j < Mb.elemNum; j++)
				Mc.elem[i * Mb.Cnum + j] = Ma.elem[i * Ma.Cnum + j] * Mb.elem[j];
		return Mc;
	}
	else{
		Matrix_t Mc(Ma.Rnum, Mb.Cnum, true);	
		for(int i = 0; i < Ma.Rnum; i++)
			Mc.elem[i] = Ma.elem[i] * Mb.elem[i];
		return Mc;
	}
}

void Matrix_t::operator*= (const Matrix_t& Mb){
	*this = *this * Mb;
}

vector<Matrix_t> Matrix_t::diagonalize(){
	assert(Rnum == Cnum);
	assert(!diag);
	vector<Matrix_t> outs;
	Matrix_t Eig(Rnum, Cnum, true);
	Matrix_t EigV(Rnum, Cnum);
	syDiag(elem, Rnum, Eig.elem, EigV.elem);
	outs.push_back(Eig);
	outs.push_back(EigV);	
	return outs;
}

vector<Matrix_t> Matrix_t::svd(){
	assert(!diag);
	vector<Matrix_t> outs;
	int min = Rnum < Cnum ? Rnum : Cnum;	//min = min(Rnum,Cnum)
	Matrix_t U(Rnum, min);
	Matrix_t S(min, min, true);
	Matrix_t VT(min, Cnum);
	myDgesvd(elem, Rnum, Cnum, U.elem, S.elem, VT.elem);
	outs.push_back(U);
	outs.push_back(S);
	outs.push_back(VT);
	return outs;
}

void Matrix_t::orthoRand(){
	assert(!diag);
	orthoRandomize(elem, Rnum, Cnum);
}

Matrix_t operator*(const Matrix_t& Ma, double a){
	Matrix_t Mb(Ma);
	vecScal(a, Mb.elem, Mb.elemNum);
	return Mb;
}

void Matrix_t::operator*= (double a){
	vecScal(a, elem, elemNum);
}

Matrix_t operator+(const Matrix_t& Ma, const Matrix_t& Mb){
	Matrix_t Mc(Ma);
	vecAdd(Mb.elem, Mc.elem, Ma.elemNum);
	return Mc;
}

void Matrix_t::operator+= (const Matrix_t& Mb){
	vecAdd(Mb.elem, elem, elemNum);
}

void Matrix_t::transpose(){
	if(!diag){
		double* oldElem = (double*)malloc(elemNum * sizeof(double));
		memcpy(oldElem, elem, elemNum * sizeof(double));
		myTranspose(oldElem, Rnum, Cnum, elem);
		free(oldElem);
	}
	int tmp = Rnum;
	Rnum = Cnum;
	Cnum = tmp;
}
