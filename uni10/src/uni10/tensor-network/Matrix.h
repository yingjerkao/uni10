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
    Matrix();

    /// @brief Create a Matrix
    ///
    /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
    /// matrix elements and set the elements to zero
    /// @param _Rnum Number of Rows
    /// @param _Cnum Number of Columns
    /// @param _diag Set \c true for diagonal matrix, defaults to \c false
    Matrix(size_t _Rnum, size_t _Cnum, bool _diag=false, bool _ongpu=false);

    /// @brief Create a Matrix
    ///
    /// Allocate memory of size <tt> Rnum * Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is \c true) for
    /// matrix elements and copy the elements from \c _elem
    Matrix(size_t _Rnum, size_t _Cnum, const double* _elem, bool _diag=false, bool src_ongpu=false);
    /// @overload
    Matrix(size_t _Rnum, size_t _Cnum, const std::vector<double>& _elem, bool _diag=false, bool src_ongpu=false);

    /// @brief Copy constructor
    Matrix(const Matrix& _m);
    /// @overload
    Matrix(const Block& _b);
#ifndef UNI10_PURE_REAL
    /// @overload
    Matrix(const CBlock& _b);
#endif

    /// @brief Destructor
    ///
    /// Destroys  Matrix and freeing all the allocated memory for matrix elements.
    ~Matrix();

    /// @brief Assign to Matrix
    ///
    /// Assigns the content of \c mat to Matrix, replacing the original content by reallocating new memory
    /// fit for \c mat.
    /// @param _m Second Matrix
    Matrix& operator=(const Matrix& _m);
    /// @overload
    Matrix& operator=(const Block& _m);

    /// @brief Access individual element
    ///
    /// Returns a reference to the element at position \c idx in Matrix. The value \c idx is a serial index
    /// counted in row-major sense from the first element (\c idx = 0) of Matrix.
    /// @note This function works similar to member function Matrix::at().
    /// @param idx Element index
    /// @return Element of Matrix at position \c idx
    double& operator[](size_t idx);

    /// @brief Access individual element
    ///
    /// Returns a reference to the element in the i-th row and j-th column of Matrix.
    /// @note The values \c i and \c j are counted from 0.
    /// @param i,j Index of Matrix
    /// @return Element at index \c (i,j) of Matrix
    double& at(size_t i, size_t j);



    double* getHostElem();



    /// @brief Copy elements
    ///
    /// Copies the first elemNum() elements from \c elem, replacing the original ones.
    /// @param elem  Matrix elements to be copied from.
    void setElem(const double* elem, bool _ongpu = false);
    /// @overload
    void setElem(const std::vector<double>& elem, bool _ongpu = false);

    /// @brief Resize Matrix
    ///
    /// Changes the number of columns and rows to the values \c row and \c col.
    /// If the size shrinks, the outsize elements are truncated. If the size grows,
    /// extra elements are filled with zero.
    /// @param row New number of rows
    /// @param col New number of columns
    Matrix& resize(size_t row, size_t col);




    /// @brief Load Matrix from a file
    ///
    /// Loads the elements of  Matrix from a binary file \c fname.
    /// @param fname Filename
    void load(const std::string& fname);

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

    bool toGPU();
    /// @brief Transpose Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& transpose();

    /// @brief Transpose a complex Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& cTranspose();

    Matrix& conj(){return *this;};

    /// @brief Multiply Matrix by a scalar and assign
    ///
    /// Performs element-wise multiplication with a scalar \c a .
    Matrix& operator*= (double a);

    /// @brief   Multiply Matrix by a second matrix and assign
    ///
    /// Performs matrix multiplication of Matrix with another Matrix \c Mb and store the results in Matrix
    Matrix& operator*= (const Block& Mb);

    /// @brief Perform addition of elements and assign
    ///
    /// Performs element by element addition and store the results in Matrix
    Matrix& operator+= (const Block & Mb);

private:
    void init(bool togpu);
    void init(const double* elem, bool _ongpu);

};
Matrix takeExp(double a, const Block& mat);
#ifndef UNI10_PURE_REAL
Matrix exp(double a, const Block& mat);
CMatrix exp(const std::complex<double>& a, const Block& mat);
#endif
Matrix exph(double a, const Block& mat);
#ifndef UNI10_PURE_REAL
Matrix exp(const Block& mat);
#endif
Matrix exph(const Block& mat);
Matrix otimes(const Block& Ma, const Block& Mb);

/// @example egM1.cpp
/// @example egM2.cpp

};  /* namespace uni10 */
#endif /* MATRIX_H */
