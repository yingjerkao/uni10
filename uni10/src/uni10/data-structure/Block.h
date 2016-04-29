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
*  @author Yun-Hsuan Chou
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
/// @brief Base class for Matrix.
///
/// A Block holds a reference to a Matrix. The Block constructor does not allocate memory. Memory allocation
/// should be done through Matrix.
///
/// If no datatype flag is specified, the constructors and member functions default to Real.
/// To explicitly declare the data type, the class provides constructors and member functions with the following syntax:
///
/// <tt> func(rflag tp, ...) </tt> for Real datatype and <tt> func(cflag tp, ...) </tt> for Complex datatype.
///
/// @see \ref Matrix, UniTensor

    class Block{
	public:

	    /*********************  OPERATOR **************************/

        /// @brief Print out Block
        ///
        /// Prints out the elements of Block
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
	    friend std::ostream& operator<< (std::ostream& os, const Block& b);
        
        /// @brief Multiplication of two matrices
        ///
        /// Performs matrix multiplication of two matrices \c Ma and \c Mb. Store the results in a new matrix
	    friend Matrix operator*(const Block& Ma, const Block& Mb); //R*R C*C R*C C*R
        
        /// @brief  Multiplication of a matrix and a scalar
        ///
        /// Performs element-wise  multiplication of \c Ma with \c  Store the results in a new matrix
	    friend Matrix operator*(Real a, const Block& Ma);
        
        /// @overload
	    friend Matrix operator*(const Block& Ma, Real a);
        /// @overload
	    friend Matrix operator*(const Complex& a, const Block& Ma);
        /// @overload
	    friend Matrix operator*(const Block& Ma, const Complex& a);
        /// @brief Addition of  two matrices
        ///
        /// Performs element by element addition of \c Ma and \c Mb. Store the results in a new matrix
	    friend Matrix operator+(const Block& Ma, const Block& Mb);
        /// @brief Compare matrices
        ///
        /// Compare all the elements in \c Ma and \c Mb, returning \c true if all the elements are the same,
        /// \c false otherwise.
        /// @return \c True if <tt> Ma == Mb </tt>, \c false otherwise.
	    friend bool operator==(const Block& m1, const Block& m2);
        /// @brief Compare matrices
        ///
        /// Compare all the elements in \c Ma and \c Mb, returning \c false if all the elements are the same,
        /// \c true otherwise.
        /// @return \c False if <tt> Ma == Mb </tt>, \c true otherwise.
	    friend bool operator!=(const Block& m1, const Block& m2){return !(m1 == m2);};
	    friend Matrix RDotR(const Block& Ma, const Block& Mb);
	    friend Matrix CDotR(const Block& Ma, const Block& Mb);
	    friend Matrix RDotC(const Block& Ma, const Block& Mb);
	    friend Matrix CDotC(const Block& Ma, const Block& Mb);
	    friend size_t lanczosEigh(Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter, Real err_tol );
	    friend size_t lanczosEigh(rflag tp, Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter, Real err_tol );
	    friend size_t lanczosEigh(cflag tp, Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter, Real err_tol );

	    /*********************  NO TYPE **************************/
        
        ///
        /// @brief Default constructor
        ///
	    Block();
        ///
        /// @brief Create a Block
        ///
        /// Create a Block of size <tt> Rnum*Cnum </tt> ( or <tt> min(Rnum, Cnum)</tt> if \c diag is\c true)
        /// without allocating memeories
        ///
        /// @param Rnum Number of rows
        /// @param Cnum Number of columns
        /// @param diag Set \c true for a diagonal matrix, defaults to \c false
        ///
        /// @overload
        Block(size_t _Rnum, size_t _Cnum, bool _diag = false);
        
        /// @brief Copy constructor
        ///
        Block(const Block& _b);
	    virtual ~Block();
        ///
        /// @brief Returns the number of rows in Block
        ///
        /// @return Number of rows
	    size_t row()const;
        ///
        /// @brief Returns the number of columns in Block
        ///
        /// @return Number of columns
	    size_t col()const;
        ///
        /// @brief Tests if Block is diagonal
        ///
        /// @return \c True if Block is diagonal, \c false otherwise
	    bool isDiag()const;
	    bool isOngpu()const;
        ///
        /// @brief Returns the number of elements of Block
        ///
        /// @return Number of elements
	    size_t elemNum()const;
        ///
        /// @brief Returns the datatype ID
        ///
        /// @return Datatype ID
	    int typeID()const;
        ///
        /// @brief Saves the Block to file
        ///
        /// @param fname Filename to save the block
	    void save(const std::string& fname)const;
        void savePrototype(const std::string& fname)const;
        ///
        /// @brief Performs QR decomposition of Block
        ///
        /// Performs QR decomposition on Block
        /// and returns a vector of two matrices \f$ [Q, R]\f$ .
        ///
        /// For \c m-by-\c n matrix \f$A\f$, it is decomposed as:
        /// \f[ A = Q\times R\f]
        ///
        /// @return A vector of matrices \f$ [Q, R]\f$.
        /// \f$Q \f$ is a unitary  matrix, and \f$R\f$ is a  upper-triangular matrix.
        ///
	    std::vector<Matrix> qr()const;
        
        ///
        /// @brief Performs RQ decomposition of Block
        ///
        ///
        /// Performs RQ decomposition on Block
        /// and returns a vector of two matrices \f$ [R, Q]\f$.
        ///
        /// @return A vector of matrices \f$ [R, Q]\f$.
        /// \f$R\f$ is a  upper-triangular matrix and \f$Q \f$ is a unitary  matrix.
	    std::vector<Matrix> rq()const;
        
        ///
        /// @brief Performs QL decomposition of Block
        ///
        ///
        /// Performs QL decomposition on Block
        /// and returns a vector of two matrices \f$ [Q, L]\f$.
        ///
        /// @return A vector of matrices \f$ [Q, L]\f$.
        /// \f$L\f$ is a  lower-triangular matrix and \f$Q \f$ is a unitary  matrix.
        std::vector<Matrix> ql()const;
        
        ///
        /// @brief Performs LQ decomposition of Block
        ///
        ///
        /// Performs LQ decomposition on Block
        /// and returns a vector of two matrices \f$ [L, Q]\f$.
        ///
        ///
        /// @return A vector of matrices \f$ [L, Q]\f$.
        /// \f$L\f$ is a  lower-triangular matrix and \f$Q \f$ is a unitary  matrix.
        std::vector<Matrix> lq()const;

        ///
        /// @brief Perform SVD decomposition 
        ///
        /// Performs singular value decomposition(SVD) on  Block
        /// and returns a vector of three matrices of SVD.
        /// For \c m-by-\c n matrix \f$A\f$, it is decomposed as:
        /// \f[ A = U \times \Sigma  \times V^\dagger \f]
        /// @return A vector of matrices \f$ [U, \Sigma, V^\dagger]\f$
        ///
        /// For an \c m by \c n Block \f$A\f$
        ///
        /// \f$U\f$ is an \c m by\c n unitary matrix (row-major).
        ///
        /// \f$\Sigma\f$ is an \c n  by \c n diagonal matrix.
        ///
        /// \f$V^\dagger \f$ is an \c n by \c m  unitary matrix (row-major).
        ///
        /// @note The operation is a wrapper of Lapack function \c Xgesvd().
        ///
	    std::vector<Matrix> svd()const;
        /// @brief  Diagonalize a General Block
        ///
        /// Diagonalizes Block and returns the eigenvalues and eigenvectors as matrices.
        ///
        /// For an \c n by \c n matrix \f$ A\f$
        ///
        /// \f[ A = U^T \times D \times U \f]
        /// \f$ D\f$ is an \c n by \c n diagonal matrix of eigenvalues.
        /// \f$ U \f$ is an \c n by \c n matrix  of right eigenvectors as row-vectors.
        /// @return A vector of matrices \f$[D, U]\f$
        /// @note Only the right eigenvectors will be given. The operation is a wrapper of Lapack function \c Xsyev().

	    std::vector<Matrix> eig()const;
        /// @brief  Diagonalize a symmetric/hermitian Block
        ///
        /// Diagonalizes Block and returns the eigenvalues and eigenvectors as matrices.
        ///
        /// For an \c n by \c n matrix \f$ A\f$
        ///
        /// \f[ A = U^T \times D \times U \f]
        /// \f$ D\f$ is an \c n by \c n diagonal matrix of eigenvalues.
        /// \f$ U \f$ is an \c n by \c n matrix  of eigenvectors as row-vectors.
        /// @return A vector of matrices \f$[D, U]\f$
        /// @note Block must be symmetric/hermitian matrix.
        /// The operation is a wrapper of Lapack function \c dsyev() for Real matrix zheev() for Complex matrix.
	    std::vector<Matrix> eigh()const;
        
        ///
        /// @brief Computes the inverse matrix of Block
        ///
        /// @return Inverse of Block
        Matrix inverse()const;
        ///
        /// @brief Computes the \f$L^2\f$-norm of the elements in Block
        ///
        /// @return The \f$L^2\f$-norm
	    Real norm()const;
        ///
        /// @brief Returns the diagonal elements of Block
        ///
        /// @return Diagonal elements in a matrix
        Matrix getDiag()const;
        ///
        ///
        /// @brief Computes the trace (sum of diagonal elements) of Block
        ///
        /// @return Trace of Block
        
	    Real trace()const;
        ///
        ///
        /// @brief Computes the sum  elements of Block
        ///
        /// @return Sum of elements of Block
        ///
	    Real sum()const;
        ///
        /// @brief Access individual element
        ///
        /// Returns the content of the element in the i-th row and j-th column of Block.
        /// Outputs the real part of the element by default
        /// To obtain the complex value, use the complex version Block::at( )const.
        /// @note The values \c i and \c j are counted from 0.
        ///
        /// @param i,j Index of Block
        /// @return Element at index \c (i,j) of Block.
	    Real at(size_t i, size_t j)const;
	    bool CelemIsNULL()const;
	    bool RelemIsNULL()const;
	    /*********************  REAL **********************/

	    Block(rflag _tp, size_t _Rnum, size_t _Cnum, bool _diag = false);
	    void save(rflag _tp, const std::string& fname)const;
	    std::vector<Matrix> qr(rflag tp)const;
	    std::vector<Matrix> rq(rflag tp)const;
	    std::vector<Matrix> ql(rflag tp)const;
	    std::vector<Matrix> lq(rflag tp)const;
	    std::vector<Matrix> svd(rflag tp)const;
	    std::vector<Matrix> eig(rflag tp)const;
	    std::vector<Matrix> eigh(rflag tp)const;
	    Matrix inverse(rflag tp)const;
	    Real norm(rflag tp)const;
	    Matrix getDiag(rflag tp)const;

	    Real trace(rflag tp)const;
	    Real sum(rflag tp)const;
	    Real operator[](size_t idx)const;
	    Real at(rflag tp, size_t i, size_t j)const;
	    Real* getElem(rflag tp = RTYPE)const;     //rename -> getRealElem() && getComplexElem();

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
	    Real norm(cflag _tp)const;
	    Matrix getDiag(cflag _tp)const;
	    Complex trace(cflag _tp)const;
	    Complex sum(cflag _tp)const;
	    Complex operator()(size_t idx)const;
	    Complex at(cflag _tp, size_t i, size_t j)const;
	    Complex* getElem(cflag _tp)const;     //rename -> getRealElem() && getComplexElem();

	    friend class UniTensor;
	    friend class Matrix;

	protected:
	    rflag r_flag;
	    cflag c_flag;
	    Real* m_elem;     // pointer to a real matrix
	    Complex* cm_elem; // pointer to a complex matrix
	    size_t Rnum;		//number of rows of the block
	    size_t Cnum;		//number of columns of the block
	    bool diag;
	    bool ongpu;
    };
    
    // Helper functions for multiplication of two blocks
    // RDotR: Real dot Real
    // CDotR: Complex dot Real
    // RDotC: Real dot Complex
    // CDotC: Complex dot Complex
    Matrix RDotR(const Block& Ma, const Block& Mb);
    Matrix CDotR(const Block& Ma, const Block& Mb);
    Matrix RDotC(const Block& Ma, const Block& Mb);
    Matrix CDotC(const Block& Ma, const Block& Mb);
    
    size_t lanczosEigh(Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter=1000, Real err_tol = 5E-15);
    size_t lanczosEigh(rflag _tp, Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter=1000, Real err_tol = 5E-15);
    size_t lanczosEigh(cflag _tp, Matrix& ori_mat, Real& E0, Matrix& psi, size_t max_iter=1000, Real err_tol = 5E-15);
};
#endif /* BLOCK_H */
