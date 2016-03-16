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
    
//! Real datatype flag
    enum rflag{
	RNULL = 0, ///< Real datatype not defined
	RTYPE = 1 ///< Real datatype defined
    };

//! Complex datatype flag
    enum cflag{
	CNULL = 0,///< Complex datatype not defined
	CTYPE = 2 ///< Complex datatype defined
    };

    class UniTensor;
    class Matrix;
/// @class Block
/// @brief Block class is the base class for Matrix.
///
/// A Block holds a reference to a Matrix. The Block constructor does not allocate memory. Memory allocation
/// should be done through Matrix.
///
/// @see \ref Matrix, UniTensor

    class Block{
	public:

	    /*********************  OPERATOR **************************/

	    friend std::ostream& operator<< (std::ostream& os, const Block& b);
	    friend Matrix operator*(const Block& Ma, const Block& Mb); //R*R C*C R*C C*R
	    friend Matrix operator*(double a, const Block& Ma);
	    friend Matrix operator*(const Block& Ma, double a);
	    friend Matrix operator*(const std::complex<double>& a, const Block& Ma);
	    friend Matrix operator*(const Block& Ma, const std::complex<double>& a);
	    friend Matrix operator+(const Block& Ma, const Block& Mb);
	    friend bool operator==(const Block& m1, const Block& m2);
	    friend bool operator!=(const Block& m1, const Block& m2){return !(m1 == m2);};

	    /*********************  NO TYPE **************************/
        
        ///
        /// @brief Default constructor
        ///
	    Block();
        ///
        /// @brief Create a Block
        ///
        /// Create a Block of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true)
        /// without allocating memeories
        ///
        /// @param _typeID Real or Complex datatype. If not present, defaults to Real.
        /// @param _Rnum Number of Rows
        /// @param _Cnum Number of Columns
        /// @param _diag Set \c true for diagonal matrix, defaults to \c false
	    Block(int _typeID, size_t _Rnum, size_t _Cnum, bool _diag = false);
        /// @overload
        Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
        
        /// @brief Copy constructor
        ///
        Block(const Block& _b);
	    virtual ~Block();
        
        /// @brief Create a Block
        ///
        /// Create a Block of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true)
        /// without allocating memeories
        ///
	    size_t row()const;
	    size_t col()const;
	    bool isDiag()const;
	    bool isOngpu()const;
	    size_t elemNum()const;
	    int typeID()const;
	    void save(const std::string& fname)const;
	    void savePrototype(const std::string& fname)const;
	    std::vector<Matrix> qr()const;
	    std::vector<Matrix> rq()const;
	    std::vector<Matrix> ql()const;
	    std::vector<Matrix> lq()const;
	    std::vector<Matrix> svd()const;
	    std::vector<Matrix> eig()const;
	    std::vector<Matrix> eigh()const;

	    Matrix inverse()const;
	    double norm()const;
	    Matrix getDiag()const;
	    std::complex<double> trace()const;
	    std::complex<double> sum()const;
	    friend size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol );
	    std::complex<double> operator[](size_t idx)const;
	    std::complex<double> at(size_t i, size_t j)const;
	    /*********************  REAL **********************/
       
	    Block(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
	    void save(rflag _tp, const std::string& fname)const;
        /// @brief Performs QR decomposition of Block
        ///
        ///
	    std::vector<Matrix> qr(rflag _tp)const;
	    std::vector<Matrix> rq(rflag _tp)const;
	    std::vector<Matrix> ql(rflag _tp)const;
	    std::vector<Matrix> lq(rflag _tp)const;
	    std::vector<Matrix> svd(rflag _tp)const;

	    std::vector<Matrix> eig(rflag _tp)const;
	    std::vector<Matrix> eigh(rflag _tp)const;
	    Matrix inverse(rflag _tp)const;
	    double norm(rflag _tp)const;
	    Matrix getDiag(rflag _tp)const;
	    double trace(rflag _tp)const;
	    double sum(rflag _tp)const;
	    double* getElem(rflag _tp)const;     //rename -> getRealElem() && getComplexElem();
	    double at(rflag _tp, size_t i, size_t j)const;
	    friend size_t lanczosEigh(rflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol );
	    /*********************  COMPLEX **********************/

	    Block(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
	    void save(cflag _tp, const std::string& fname)const;
	    std::vector<Matrix> qr(cflag _tp)const;
	    std::vector<Matrix> rq(cflag _tp)const;
	    std::vector<Matrix> ql(cflag _tp)const;
	    std::vector<Matrix> lq(cflag _tp)const;
	    std::vector<Matrix> svd(cflag _tp)const;

	    std::vector<Matrix> eig(cflag _tp)const;
	    std::vector<Matrix> eigh(cflag _tp)const;
	    Matrix inverse(cflag _tp)const;
	    double norm(cflag _tp)const;
	    Matrix getDiag(cflag _tp)const;
	    std::complex<double> trace(cflag _tp)const;
	    std::complex<double> sum(cflag _tp)const;
	    std::complex<double>* getElem(cflag _tp)const;     //rename -> getRealElem() && getComplexElem();
	    std::complex<double> at(cflag _tp, size_t i, size_t j)const;
	    friend size_t lanczosEigh(cflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter, double err_tol );
	    void printElemIsNULL() const;
	    /******************Friend funcs*******************/

	    friend void RtoC(Block& mat);
	    friend void RtoC(UniTensor& UniT);
	    friend Matrix exph(double a, const Block& mat);

	    /*****************************************************/

	    double* getElem()const;     //rename -> getRealElem() && getComplexElem();
	    //friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);

	    /**********************************************************/
	    friend class UniTensor;
	    friend class Matrix;
	
	    /********************************************************************************/

	protected:
	    rflag r_flag;
	    cflag c_flag;
	    double* m_elem;     // pointer to a real matrix
	    std::complex<double>* cm_elem; // pointer to a complex matrix
	    size_t Rnum;		//number of rows of the block
	    size_t Cnum;		//number of columns of the block
	    bool diag;
	    bool ongpu;
    };
    void RtoC(Block& mat);
    // size_t lanczosEigh(Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter=1000, double err_tol = 5E-15);
    // size_t lanczosEigh(rflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter=1000, double err_tol = 5E-15);
    // size_t lanczosEigh(cflag _tp, Matrix& ori_mat, double& E0, Matrix& psi, size_t max_iter=1000, double err_tol = 5E-15);
};
#endif /* BLOCK_H */
