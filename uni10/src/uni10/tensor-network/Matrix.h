/****************************************************************************
*  @file Matrix.h
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
#include <cstdio>
#include <uni10/data-structure/Block.h>
#include <uni10/data-structure/CBlock.h>
//Type of Matrix
namespace uni10{

/// @class Matrix
/// @brief The Matrix class defines a common matrix
///
/// Matrix is an auxilliary class used to extract block information on UniTensor and perform linear algebra
/// operations. A symmetric tensor contains Block's with corresponding Qnum's. Each block is a Matrix with
/// tensor elements.
///
/// The Matrix follows the C convention that the memory storage is row-major and indices start from 0.
class Matrix:public Block {
public:
    ///@brief Default constructor
    ///
	  /*********************  OPERATOR **************************/	    

    Matrix& operator=(const Matrix& _m);
    Matrix& operator=(const Block& _m);
    Matrix& operator*= (double a);
    Matrix& operator*= (std::complex<double> a);
    Matrix& operator*= (const Block& Mb);
    Matrix& operator+= (const Block& Mb);
    std::complex<double> operator[](size_t idx); //&

    /********************* going move **************************/	    

    Matrix(muType tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    Matrix(const CBlock& _b);
    void assign(muType _tp, size_t _Rnum, size_t _Cnum);

    /*********************  NO TYPE **************************/	    

    Matrix();
    Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    Matrix(const Matrix& _m);
    Matrix(const Block& _b);
    Matrix(const std::string& fname);
    ~Matrix();
    void identity();
    void set_zero();
    void randomize();
    void orthoRand();
    Matrix& transpose();
    Matrix& cTranspose();
    Matrix& conj();
    Matrix& resize(size_t row, size_t col);
    double max(bool _ongpu=false);
    void load(const std::string& fname);
    void assign(size_t _Rnum, size_t _Cnum);
    bool toGPU();

    /*********************  REAL **********************/

    Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag=false, bool src_ongpu=false);
    Matrix(size_t _Rnum, size_t _Cnum, const std::complex<double>* _elem, bool _diag=false, bool src_ongpu=false);
    Matrix(rflag _tp, const std::string& fname);
    Matrix(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    void setElem(const double* elem, bool _ongpu = false);
    void setElem(const std::vector<double>& elem, bool _ongpu = false);
    void identity(rflag _tp);
    void set_zero(rflag _tp);
    void randomize(rflag _tp);
    void orthoRand(rflag _tp);
    Matrix& transpose(rflag _tp);
    Matrix& cTranspose(rflag _tp);
    Matrix& conj(rflag _tp);
    Matrix& resize(rflag _tp, size_t row, size_t col);
    double max(rflag _tp, bool _ongpu=false);
    void assign(rflag _tp, size_t _Rnum, size_t _Cnum);
    bool toGPU(rflag _tp);
    double& at(rflag _tp, size_t i); //&
    double* getHostElem(rflag _tp);

    /*********************  COMPLEX **********************/
    
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag=false, bool src_ongpu=false);
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector< std::complex<double> >& _elem, bool _diag=false, bool src_ongpu=false);
    Matrix(cflag _tp, const std::string& fname);
    Matrix(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    void setElem(const std::complex<double>* elem, bool _ongpu = false);
    void setElem(const std::vector< std::complex<double> >& elem, bool _ongpu = false);
    void identity(cflag _tp);
    void set_zero(cflag _tp);
    void randomize(cflag _tp);
    void orthoRand(cflag _tp);
    Matrix& transpose(cflag _tp);
    Matrix& cTranspose(cflag _tp);
    Matrix& conj(cflag _tp);
    Matrix& resize(cflag _tp, size_t row, size_t col);
    double max(cflag _tp, bool _ongpu=false);
    void assign(cflag _tp, size_t _Rnum, size_t _Cnum);
    bool toGPU(cflag _tp);
    std::complex<double>& at(cflag _tp, size_t i); //&
    std::complex<double>* getHostElem(cflag _tp);
   
    /*****************************************************/

    //delete
    double* getHostElem();
    std::complex<double> at(size_t i, size_t j); //&

/********************************************************************************/
private:
    /********************* going move **************************/	    
    void init(bool togpu, muType tp);
    /*********************  NO TYPE **********************/
    void melemFree();
    void setMelemBNULL();
    void init(const double* _m_elem, const std::complex<double>* _cm_elem, bool src_ongpu);
    /*********************  REAL **********************/
    void init(const double* elem, bool _ongpu);
    void init(rflag _tp, bool togpu);
    /*********************  COMPLEX **********************/
    void init(const std::complex<double>* elem, bool _ongpu);
    void init(cflag _tp, bool togpu);
    /***********************************8*****************/
};

Matrix takeExp(double a, const Block& mat);
Matrix exp(double a, const Block& mat);
Matrix exp(const std::complex<double>& a, const Block& mat);
Matrix exph(double a, const Block& mat);
Matrix exp(const Block& mat);
Matrix exph(const Block& mat);
Matrix otimes(const Block& Ma, const Block& Mb);

/// @example egM1.cpp
/// @example egM2.cpp

};  /* namespace uni10 */
#endif /* MATRIX_H */
