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

    friend void RtoC(Matrix& mat);
    friend void RAddR(Matrix& Ma, const Matrix& Mb);
    friend void CAddR(Matrix& Ma, const Matrix& Mb);
    friend void RAddC(Matrix& Ma, const Matrix& Mb);
    friend void CAddC(Matrix& Ma, const Matrix& Mb);

    /*********************  developping  **********************/
    
    Real absMax(bool _ongpu=false);
    Real absMax(rflag tp, bool _ongpu=false);

    Matrix& maxNorm();
    Matrix& maxNorm(rflag tp);

    Matrix& absMaxNorm();
    Matrix& absMaxNorm(rflag tp);

    Real* getHostElem();
    Real* getHostElem(rflag _tp);
    Complex* getHostElem(cflag _tp);

    bool toGPU();
    
	  /*********************  OPERATOR **************************/
    /// @brief Assigns to Matrix
    ///
    /// Assigns the content of \c mat to Matrix, replacing the original content by reallocating new memory
    /// fit for \c mat.
    /// @param _m Second Matrix
    Matrix& operator=(const Matrix& _m);
    /// @overload
    Matrix& operator=(const Block& _m);

    /// @brief Multiply Matrix by a scalar and assign
    ///
    /// Performs element-wise multiplication with a Real scalar \c a .
    Matrix& operator*= (Real a);

    /// @brief Multiply Matrix by a scalar and assign
    ///
    /// Performs element-wise multiplication with a Complex scalar \c a .
    Matrix& operator*= (Complex a);

    /// @brief   Multiply Matrix by a second matrix and assign
    ///
    /// Performs matrix multiplication of Matrix with another Matrix \c Mb and store the results in Matrix
    Matrix& operator*= (const Block& Mb);

    /// @brief Perform addition of elements and assign
    ///
    /// Performs element by element addition and store the results in Matrix

    Matrix& operator+= (const Block& Mb);

    /*********************  NO TYPE **************************/
    ///@brief Default constructor
    ///
    Matrix();

    /// @brief Create a Matrix
    ///
    /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
    /// matrix elements and set the elements to zero
    /// @param _Rnum Number of Rows
    /// @param _Cnum Number of Columns
    /// @param _diag Set \c true for diagonal matrix, defaults to \c false
    Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    Matrix(std::string tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);
    /// @brief Copy constructor
    Matrix(const Matrix& _m);
    /// @overload
    Matrix(const Block& _b);
    // #### new
    Matrix(const std::string& fname);
    /// @brief Destructor
    ///
    /// Destroys Matrix and freeing all the allocated memory for matrix elements.
    ~Matrix();

    /// @brief Set to identity
    ///
    /// Sets the diagonal matrix elements to 1.
    void identity();

    /// @brief Assign zero elements
    ///
    /// Assigns zeros to the elements of Matrix.
    void set_zero();

    /// @brief Assign random elements
    ///
    /// Assigns random values between [0, 1) to  elements of Matrix.
    void randomize();

    /// @brief Assign elements
    ///
    /// Assigns randomly generated orthogonal bases to the elements of  Matrix.
    /// \c Nr = row() and \c Nc = col().
    /// If the <tt> Nr < Nc </tt>, randomly generates \c Nr orthogonal basis row vectors of dimension \c Nc.
    /// If the <tt> Nr > Nc </tt>, randomly generates \c Nc orthogonal basis column vectors of dimension \c Nr.
    void orthoRand();

    Matrix& normalize();

    /// @brief Transpose Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& transpose();

    /// @brief Transpose a complex Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& cTranspose();
    Matrix& conj();

    /// @brief Resize Matrix
    ///
    /// Changes the number of columns and rows to the values \c row and \c col.
    /// If the size shrinks, the outsize elements are truncated. If the size grows,
    /// extra elements are filled with zero.
    /// @param row New number of rows
    /// @param col New number of columns
    Matrix& resize(size_t row, size_t col);

    /// @brief Returns the maximum element in Matrix
    ///
    /// Returns the maximum element in Matrix
    Real max(bool _ongpu=false);

    void assign(size_t _Rnum, size_t _Cnum);

    /// @brief Load Matrix from a file
    ///
    /// Loads the elements of  Matrix from a binary file \c fname.
    /// @param fname Filename
    void load(const std::string& fname);

    /*********************  REAL **********************/

    /// @brief Create a Matrix
    ///
    /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
    /// matrix elements and copy the elements from \c _elem
    Matrix(size_t _Rnum, size_t _Cnum, const Real* _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    /// @overload
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector<Real>& _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    // @overload
    Matrix(rflag _tp, const std::string& fname);
    // @overload
    Matrix(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);

    /// @brief Copy elements
    ///
    /// Copies the first elemNum() elements from \c elem, replacing the original ones.
    /// @param elem  Matrix elements to be copied from.
    void setElem(const Real* elem, bool src_ongpu = false);

    /// @overload
    void setElem(const std::vector<Real>& elem, bool src_ongpu = false);

    /// @brief Set to identity
    ///
    /// Sets the diagonal matrix elements to 1.
    void identity(rflag _tp);

    /// @brief Assign zero elements
    ///
    /// Assigns zeros to the elements of Matrix.
    void set_zero(rflag _tp);

    /// @brief Assign random elements
    ///
    /// Assigns random values between [0, 1) to  elements of Matrix.
    void randomize(rflag _tp);

    /// @brief Assign elements
    ///
    /// Assigns randomly generated orthogonal bases to the elements of  Matrix.
    /// \c Nr = row() and \c Nc = col().
    /// If the <tt> Nr < Nc </tt>, randomly generates \c Nr orthogonal basis row vectors of dimension \c Nc.
    /// If the <tt> Nr > Nc </tt>, randomly generates \c Nc orthogonal basis column vectors of dimension \c Nr.
    void orthoRand(rflag _tp);

    Matrix& normalize(rflag tp);

    /// @brief Transpose Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& transpose(rflag _tp);
    

    /// @brief Transpose a complex Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& cTranspose(rflag _tp);
    Matrix& conj(rflag _tp);

    /// @brief Resize Matrix
    ///
    /// Changes the number of columns and rows to the values \c row and \c col.
    /// If the size shrinks, the outsize elements are truncated. If the size grows,
    /// extra elements are filled with zero.
    /// @param row New number of rows
    /// @param col New number of columns
    Matrix& resize(rflag _tp, size_t row, size_t col);

    /// @brief Returns the maximum element in Matrix
    ///
    /// Returns the maximum element in Matrix
    Real max(rflag _tp, bool _ongpu=false);

    /// @brief Access individual element
    ///
    /// Returns a reference to the element in the i-th row and j-th column of Matrix.
    /// @note The values \c i and \c j are counted from 0.
    /// @param i,j Index of Matrix
    /// @return Element at index \c (i,j) of Matrix
    Real& at(size_t i, size_t j); //&

    /// @brief Access individual element
    ///
    /// Returns a complex value to the element at position \c idx in Matrix. The value \c idx is a serial index
    /// counted in row-major sense from the first element (\c idx = 0) of Matrix.
    /// @note This function works similar to member function Matrix::at().
    /// @param idx Element index
    /// @return Element of Matrix at position \c idx
    Real& at(rflag _tp, size_t i); //&
    Real& operator[](size_t idx); //&

    void assign(rflag _tp, size_t _Rnum, size_t _Cnum);

    bool toGPU(rflag _tp);

    /*********************  COMPLEX **********************/

    /// @brief Create a Matrix
    ///
    /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
    /// matrix elements and copy the elements from \c _elem
    Matrix(size_t _Rnum, size_t _Cnum, const Complex* _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);
    /// @overload
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector< Complex >& _elem, bool _diag=false, bool _ongpu=false, bool src_ongpu=false);

    Matrix(cflag _tp, const std::string& fname);

    Matrix(cflag _tp, size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);

    /// @brief Copy elements
    ///
    /// Copies the first elemNum() elements from \c elem, replacing the original ones.
    /// @param elem  Matrix elements to be copied from.
    void setElem(const Complex* elem, bool src_ongpu = false);

    /// @overload
    void setElem(const std::vector< Complex >& elem, bool src_ongpu = false);

    /// @brief Set to identity
    ///
    /// Sets the diagonal matrix elements to 1.
    void identity(cflag _tp);

    /// @brief Assign zero elements
    ///
    /// Assigns zeros to the elements of Matrix.
    void set_zero(cflag _tp);

    /// @brief Assign random elements
    ///
    /// Assigns random values between [0, 1) to  elements of Matrix.
    void randomize(cflag _tp);

    /// @brief Assign elements
    ///
    /// Assigns randomly generated orthogonal bases to the elements of  Matrix.
    /// \c Nr = row() and \c Nc = col().
    /// If the <tt> Nr < Nc </tt>, randomly generates \c Nr orthogonal basis row vectors of dimension \c Nc.
    /// If the <tt> Nr > Nc </tt>, randomly generates \c Nc orthogonal basis column vectors of dimension \c Nr.
    void orthoRand(cflag _tp);

    Matrix& normalize(cflag tp);

    /// @brief Transpose Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& transpose(cflag _tp);

    /// @brief Transpose a complex Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& cTranspose(cflag _tp);
    Matrix& conj(cflag _tp);

    /// @brief Resize Matrix
    ///
    /// Changes the number of columns and rows to the values \c row and \c col.
    /// If the size shrinks, the outsize elements are truncated. If the size grows,
    /// extra elements are filled with zero.
    /// @param row New number of rows
    /// @param col New number of columns
    Matrix& resize(cflag _tp, size_t row, size_t col);

    /// @brief Access individual element
    ///
    /// Returns a complex value to the element at position \c idx in Matrix. The value \c idx is a serial index
    /// counted in row-major sense from the first element (\c idx = 0) of Matrix.
    /// @note This function works similar to member function Matrix::at().
    /// @param idx Element index
    /// @return Element of Matrix at position \c idx
    Complex& operator()(size_t idx); //&
    Complex& at(cflag _tp, size_t i); //&

    void assign(cflag _tp, size_t _Rnum, size_t _Cnum);

    bool toGPU(cflag _tp);

    
    /*****************************************************/

private:
    /*********************  NO TYPE **********************/

    void MelemFree();
    void setMelemBNULL();
    void init(const Real* _m_elem, const Complex* _cm_elem, bool src_ongpu);

    /*********************  REAL **********************/

    void init(const Real* elem, bool _ongpu);
    void init(rflag tp, bool togpu);

    /*********************  COMPLEX **********************/

    void init(const Complex* elem, bool _ongpu);
    void init(cflag tp, bool togpu);

    /***********************************8*****************/
};

void RtoC(Matrix& mat);
void RAddR(Matrix& Ma, const Matrix& Mb);
void CAddR(Matrix& Ma, const Matrix& Mb);
void RAddC(Matrix& Ma, const Matrix& Mb);
void CAddC(Matrix& Ma, const Matrix& Mb);
Matrix takeExp(Real a, const Block& mat);
Matrix exph(Real a, const Block& mat);
Matrix exph(rflag tp, Real a, const Block& mat);
Matrix exph(cflag tp,Real a, const Block& mat);
Matrix exp(Real a, const Block& mat);
Matrix exp(const Complex& a, const Block& mat);
Matrix exp(const Block& mat);
Matrix exph(const Block& mat);

Matrix otimes(const Block& Ma, const Block& Mb);


/// @example egM1.cpp
/// @example egM2.cpp

};  /* namespace uni10 */
#endif /* MATRIX_H */
