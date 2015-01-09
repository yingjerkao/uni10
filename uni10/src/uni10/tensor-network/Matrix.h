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
//Type of Matrix
//#include <uni10/tensor-network/UniTensor.h>
namespace uni10 {

/// @class Matrix
/// @brief The Matrix class defines a common matrix
///
/// Matrix is an auxilliary class used to extract block information on UniTensor and perform linear algebra
/// operations. A symmetric tensor contains Block's with corresponding Qnum's. Each block is a Matrix with
/// tensor elements.
///
/// The Matrix follows the C convention that the memory storage is row-major and indices start from 0.
class Matrix {
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
    
    /// @brief Access the number of rows
    ///
    /// Returns the number of rows in Matrix
    /// @return number of rows
    size_t row()const;
    
    /// @brief Access the number of columns
    ///
    /// Returns the number of colums in Matrix
    /// @return number of columns
    size_t col()const;
    
    /// @brief Test whether Matrix is diagonal
    ///
    /// Returns whether  Matrix is diagonal
    /// @return \c True if Matrix is diagonal; \c false otherwise
    bool isDiag()const {
        return diag;
    };
    
    bool isOngpu()const {
        return ongpu;
    };
    
    /// @brief Access the number of elements
    ///
    /// Returns the number of allocated elements. If diagonal, the return value is equal to the number of
    /// diagonal elements.
    /// @return  Number of elements in Matrix
    size_t elemNum()const;
    
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
    
    /// @brief Access elements
    ///
    /// Returns a pointer to the first element of Matrix.
    /// @return Pointer to Matrix elements
    double* getElem()const;
    
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
    
    /// @brief Output Matrix to a file
    ///
    /// Writes the elements of  Matrix to a binary file \c fname.
    /// @param fname Filename
    void save(const std::string& fname);
    
    /// @brief Output Matrix to a file
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
    
    /// @brief Transpose Matrix
    ///
    /// Transposes the elements of the Matrix. Exchange the row and column numbers.
    Matrix& transpose();
    
    /// @brief  Diagonalize Matrix
    ///
    /// Diagonalizes  Matrix and returns the eigenvalues and eigenvectors as matrices.
    ///
    /// For an \c n by \c n matrix \f$ A\f$
    ///
    /// \f[ A = U^T \times D \times U \f]
    /// \f$ D\f$ is an \c n by \c n diagonal matrix of eigenvalues.
    /// \f$ U \f$ is an \c n by \c n matrix  of eigenvectors as row-vectors.
    /// @return A vector of Matrices \f$[D, U]\f$
    /// @note Matrix must be symmetric matrix. The operation is a wrapper of Lapack function \c dsyev().
    std::vector<Matrix> eigh()const;
    
    /// @brief Perform SVD
    ///
    /// Performs singular value decomposition(SVD) on  Matrix and returns a vector of three matrices of SVD.
    /// For \c m by \c n matrix \f$A\f$, it is decomposed as:
    /// \f[ A = U \times \Sigma  \times VT \f]
    
    /// @return A vector of Matrices \f$ [U, \Sigma, VT]\f$
    /// For an \c m by \c n Matrix \f$A\f$
    ///
    /// \f$U\f$ is \c m by\c n row-major matrix.
    ///
    /// \f$\Sigma\f$ is an \c n  by \c n diagonal matrix.
    ///
    /// \f$VT\f$ is an \c n by \c m row-major matrix.
    ///
    /// @note The operation is a wrapper of Lapack function \c dgesvd().
    std::vector<Matrix> svd()const;
    
    /// @brief Find the lowest eigenvalue and the  eigenvector
    ///
    /// Find the lowest eigenvalue and the  eigenvector of Matrix by Lanczos
    /// @param E0 Lowest eigenvalue on success.
    /// @param psi initial vector (diagonal Matrix) for Lanczos iteration process, and output the final
    /// eigenvector on success
    /// @param  max_iter maximum number of iterations
    /// @param err_tol tolerance of error for the iteration process
    /// @return number of iterations
    size_t lanczosEigh(double& E0, Matrix& psi, size_t max_iter=200, double err_tol = 5E-15);
    
    /// @brief Trace of  Matrix
    ///
    /// Performs trace of  Matrix.
    /// @return Trace  of  Matrix.
    double trace();
    
    /// @brief \f$L^2\f$ norm
    ///
    /// Computes the \f$L^2\f$ norm of Matrix
    /// @return  norm of Matrix
    double norm();
    
    /// @brief Sum of the elements
    ///
    /// Sums over the elements of the Matrix
    /// @return Sum
    double sum();
    
    /// @brief Multiply Matrix by a scalar and assign
    ///
    /// Performs element-wise multiplication with a scalar \c a .
    Matrix& operator*= (double a);
    
    /// @brief   Multiply Matrix by a second matrix and assign
    ///
    /// Performs matrix multiplication of Matrix with another Matrix \c Mb and store the results in Matrix
    Matrix& operator*= (const Matrix& Mb);
    /// @brief Perform addition of elements and assign
    ///
    /// Performs element by element addition and store the results in Matrix
    Matrix& operator+= (const Matrix& Mb);
    
    /// @brief Exponetiate a matrix
    ///
    /// Exponetiate matrix \c mat as \f$exp(a \times mat)\f$.
    /// @param a exponent parameter
    /// @param mat input matrix
    friend Matrix takeExp(double a, const Matrix& mat);
    
    /// @brief Multiplication of two matrices
    ///
    /// Performs matrix multiplication of two matrices \c Ma and \c Mb. Store the results in a new matrix
    friend Matrix operator* (const Matrix& Ma, const Matrix& Mb);
    
    /// @brief  Multiplication of a matrix and a scalar
    ///
    /// Performs element-wise  multiplication of \c Ma with \c  Store the results in a new matrix
    friend Matrix operator*(const Matrix& Ma, double a);
    /// @overload
    friend Matrix operator*(double a, const Matrix& Ma) {
        return Ma * a;
    };
    /// @brief Perform addition of elements
    ///
    /// Performs element by element addition of \c Ma and \c Mb. Store the results in a new matrix
    friend Matrix operator+(const Matrix& Ma, const Matrix& Mb);
    
    /// @brief Compare Matrices
    ///
    /// Compare all the elements in \c Ma and \c Mb, returning \c true if all the elements are the same,
    /// \c false otherwise.
    /// @return \c True if <tt> Ma == Mb </tt>, \c false otherwise.
    friend bool operator== (const Matrix& Ma, const Matrix& Mb);
    /// @brief Print out Matrix
    ///
    /// Prints out the elements of Matrix
    ///
    /// For a 2 x 3 matrix \c M,
    /// \code
    /// 2 x 3 = 6
    ///
    /// -0.254 -0.858 -0.447
    ///
    /// 0.392  0.331 -0.859
    /// \endcode
    /// The number of elements is 2 x 3 = 6 and followed by a 2 by 3 matrix of elements.
    friend std::ostream& operator<< (std::ostream& os, const Matrix& b);
    bool toGPU();
private:
    void init(bool togpu);
    void init(const double* elem, bool _ongpu);
    size_t Rnum;        //number of rows of the block
    size_t Cnum;        //number of columns of the block
    double* m_elem;
    size_t m_elemNum;
    bool diag;
    bool ongpu;
};
Matrix takeExp(double a, const Matrix& mat);
Matrix otimes(const Matrix& Ma, const Matrix& Mb);
    
/// @example egM1.cpp
/// @example egM2.cpp

};  /* namespace uni10 */
#endif /* MATRIX_H */
