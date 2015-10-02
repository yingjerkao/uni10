/****************************************************************************
*  @file Block.h
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

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
*  @brief Header file for Block class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef BLOCK_H
#define BLOCK_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cstdint>
#include <cstdio>
#include <vector>
#include <uni10/datatype.hpp>
#include <stdexcept>
#include <complex>

namespace uni10{
	    
    enum muType {
	EMPTY = -1, ///<Defines an incoming Bond
	RL = 0,  ///<Defines an outgoing Bond
	CX = 1
    };
    
    class UniTensor;
    class Matrix;
    class CMatrix;
    class Block{
	public:
	    friend std::ostream& operator<< (std::ostream& os, const muType& tp);

	    /********************* verified **************************/	    

	    Block();
	    Block(const Block& _b);
	    Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
	    Block(muType _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
	    virtual ~Block();
	    size_t row()const;
	    size_t col()const;
	    bool isDiag()const;
	    bool isOngpu()const;
	    size_t elemNum()const;
	    void save(const std::string& fname)const;
	    double norm()const;
	    muType getType()const;
	    std::vector<Matrix> qr()const;
	    std::vector<Matrix> rq()const;
	    std::vector<Matrix> ql()const;
	    std::vector<Matrix> lq()const;
	    std::vector<Matrix> svd()const;
	    Matrix inverse()const;
	    double* getElem()const;     //rename -> getRealElem() && getComplexElem();
	    double* getRealElem()const;
	    std::complex<double>* getComplexElem()const;
	    friend void RtoC(Block& mat);
	    friend void RtoC(UniTensor& UniT);
	    friend Matrix operator*(const Block& Ma, const Block& Mb); //R*R C*C R*C C*R
	    friend Matrix operator*(double a, const Block& Ma);
	    friend Matrix operator*(const Block& Ma, double a);
	    friend Matrix operator*(const std::complex<double>& a, const Block& Ma);
	    friend Matrix operator*(const Block& Ma, const std::complex<double>& a);
	    friend Matrix operator+(const Block& Ma, const Block& Mb);
	    friend bool operator==(const Block& m1, const Block& m2);
	    friend bool operator!=(const Block& m1, const Block& m2){return !(m1 == m2);};
	    std::complex<double> trace()const;
	    std::complex<double> sum()const;
	    std::complex<double> operator[](size_t idx)const;
	    std::complex<double> at(size_t i, size_t j)const;
	    std::vector<Matrix> eig()const;
	    std::vector<Matrix> eigh()const;
	    friend size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol );
	    friend std::ostream& operator<< (std::ostream& os, const Block& b);
	    friend Matrix exph(double a, const Block& mat);
	    Matrix getDiag()const;
	    
	    /**********************************************************/	    
	    
	    
	    friend class UniTensor;
	    friend class CUniTensor;
	    friend class CBlock;
	    friend class Matrix;
	    friend class CMatrix;
	    /********************************************************************************/	    
	    //friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);

	protected:
	    muType m_type;
	    double* m_elem;
	    std::complex<double>* cm_elem;
	    size_t Rnum;		//number of rows of the block
	    size_t Cnum;		//number of columns of the block
	    bool diag;
	    bool ongpu;
    };
    void RtoC(Block& mat);
    size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter=1000, double err_tol = 5E-15);
};
#endif /* BLOCK_H */
