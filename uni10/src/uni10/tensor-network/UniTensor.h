/****************************************************************************
 *  @file CMakeLists.txt
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
 *  @brief Header file for UniTensor class
 *  @author Yun-Da Hsieh
 *  @date 2014-05-06
 *  @since 0.1.0
 *
 *****************************************************************************/
#ifndef UNITENSOR_H
#define UNITENSOR_H
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <assert.h>
#include <stdint.h>
#include <sstream>
#include <stdexcept>
#define DOUBLE double

#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/data-structure/Block.h>
#include <uni10/tensor-network/Matrix.h>
/// @brief Uni10 - the Universal Tensor %Network Library


namespace uni10 {
    
    ///@class UniTensor
    ///@brief The UniTensor class defines the symmetric tensors
    ///
    /// A UniTensor consists of Bond's carrying quantum numbers Qnum's. The tensor elements are organized as
    /// quantum number blocks. The Qnum's on the Bonds defines the size of the Qnum blocks and the rank of
    /// UniTensor is defined by the number of Bond's.\par
    /// Each Bond carries a label. Labels are used to manipulate tensors, such as in permute, partialTrace and
    /// contraction. \par
    /// Operations on tensor elements is pefromed through  getBlock and putBlock functions to take out/put in
    /// block elements out as a Matrix.
    /// @see Qnum, Bond, Matrix
    /// @example egQ1.cpp

    class UniTensor {
    public:
        
        ///
        /// @brief Default Constructor
        ///
        UniTensor();
        /// @brief Create a rank-0 UniTensor with value \c val
        ///
        /// @param val Value of the scalar
        UniTensor(double val);
        
        /// @brief Create a UniTensor from a file
        ///
        /// @param fname Filename to be read in
        UniTensor(const std::string& fname);
        
        /// @brief Create a UniTensor from a list of Bond's
        /// @param _bonds List of bonds
        /// @param _name Name of the tensor, defaults to ""
        UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
        
        /// @brief Create a UniTensor from a list of Bond's and assign labels
        /// @param _bonds List of bonds
        /// @param labels Labels for \c _bonds
        /// @param _name Name of the tensor, defaults to ""
        UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        
        /// @brief Create a UniTensor from a list of Bond's and assign labels
        /// @param _bonds List of bonds
        /// @param labels Labels for \c _bonds
        /// @param _name Name of the tensor, defaults to ""
        UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        
        /// @brief Copy constructor
        UniTensor(const UniTensor& UniT);
        /// @brief Destructor
        ///
        /// Destroys the UniTensor and freeing all allocated memory
        ~UniTensor();
    

        /// @brief Copy content
        ///
        /// Assigns new content to the UniTensor from \c UniT, replacing the original contents
        /// @param UniT Tensor to be copied
        ///
        UniTensor& operator=(const UniTensor& UniT);
        
        /// @brief Assign bonds
        ///
        /// Reconstructs the tensor with given bonds, replacing the original ones and clear the content of
        /// UniTensor
        /// @param _bond array of bonds
        UniTensor& assign(const std::vector<Bond>& _bond);
        
        /// @brief Assign labels to bonds in UniTensor
        ///
        /// Assigns the labels \c newLabels to each bond of  UniTensor, replacing the origin labels on the bonds.
        /// @param newLabels array of labels
        void setLabel(const std::vector<int>& newLabels);
        
        /// @overload
        void setLabel(int* newLabels);
        
        /// @brief Access labels
        ///
        /// Returns the labels of the bonds in UniTensor
        ///
        /// @return List of labels
        std::vector<int> label()const;
        
        /// @brief Access label
        ///
        /// Access the label of Bond \c idx
        /// @param idx Bond index
        /// @return Label of Bond \c idx
        int label(size_t idx)const;
        
        /// @brief Access the number of bonds
        ///
        /// @return Number of bonds
        size_t bondNum()const;
        
        /// @brief Access the number of incoming bonds
        ///
        /// @return Number of incoming bonds
        size_t inBondNum()const;
        
        /// @brief Access bonds
        ///
        /// Returns the bonds in UniTensor
        /// @return List of bonds
        std::vector<Bond> bond()const;
        
        /// @brief Access bond
        ///
        /// Returns the bond at the position \c idx
        /// @param idx Position of the bond being retrieved
        /// @return A bond
        Bond bond(size_t idx)const;
        
        /// @brief Access the number of elements
        ///
        /// Returns the number of total elements of the blocks.
        /// @return  Number of elements
        size_t elemNum()const;
        
        /// @brief Access raw elements
        ///
        /// Returns the elements of UniTensor in the non-block-diagonal form (raw elements) as a Matrix.
        /// The row(or column) bases of the elements are defined by the incoming bonds(or outgoing) bonds.
        /// @return Matrix of raw elements
        Matrix getRawElem()const;
        
        /// @brief Assign raw elements
        ///
        /// Assigns raw elements in \c rawElem to UniTensor
        /// The raw elements are organized in the basis given in the bonds.
        ///
        /// This function will reorganize the raw elements into the block-diagonal form.
        /// @param rawElem array of raw elements
        void setRawElem(const std::vector<double>& rawElem);
        /// @overload
        void setRawElem(const double* rawElem);
        
        /// @brief Access single element
        ///
        /// Returns the element at position specified by the indices \c idxs.
        /// The size of the array \c idxs is equal to the total number of bonds
        /// @param idxs  An array of indices
        /// @return The element at indices \c idxs.
        double at(const std::vector<int>& idxs)const;
        /// @overload
        double at(const std::vector<size_t>& idxs)const;
        
        /// @brief Access the number of blocks
        ///
        /// Returns the number of blocks
        /// @return The number of blocks
        size_t blockNum()const;
        
        /// @brief Access block quantum numbers
        ///
        /// Returns the quantum numbers for all blocks in UniTensor.
        /// The return array of quantum numbers is in the ascending order defined in Qnum.
        /// @return Array of Qnum's
        std::vector<Qnum> blockQnum()const;
        /// @brief Access block quantum number
        ///
        /// Returns the quantum number for block \c idx in UniTensor.
        /// Blocks are orderd in the ascending order of Qnum
        /// @param idx Block index
        /// @return Quantum number of block \c idx
        Qnum blockQnum(size_t idx)const;
        
        /// @brief Access block elements
        ///
        /// Returns the tensor element blocks of Qnums as
        ///  a map from a composite Qnum to a corresponding element block as Matrix
        ///
        /// @return   Map from Qnum to Matrix
        std::map<Qnum, Matrix> getBlocks()const;
        
        /// @brief Get elements in a block
        ///
        /// Returns a Matrix of Qnum(0) block elements. If the \c diag flag is set,
        /// only  diagonal elements in the block will be copied to a diagonal Matrix.
        /// @param diag Set \c true to save only the diagonal elements
        /// @return A Matrix of Qnum(0) block
        Matrix getBlock(bool diag = false)const;
        
        /// @brief Get elements in a block
        ///
        /// Returns the block elements of a quantum number \c qnum as a Matrix. If the \c diag flag is set,
        /// only  diagonal elements in the block will be copied to a diagonal Matrix.
        /// @param diag Set \c true to save only the diagonal elements
        /// @return A Matrix of \c qnum block
        Matrix getBlock(const Qnum& qnum, bool diag = false)const;
        
        /// @brief Assign elements to a block
        ///
        /// Assigns elements of the  matrix \c mat to the  block of quantum number \c qnum, replacing the origin
        /// elements. \par
        /// If \c mat is diagonal,  all the off-diagonal elements are set to zero.
        /// @param mat The matrix elements to be assigned
        void putBlock(const Matrix& mat);
            
        /// @brief Assign elements to a block
        ///
        /// Assigns elements of the  matrix \c mat to the  Qnum(0) block, for non-symmetry tensors.
        ///
        /// If \c mat is diagonal,  all the off-diagonal elements are set to zero.
        /// @param qnum quantum number of the block
        /// @param mat The matrix elements to be assigned
        void putBlock(const Qnum& qnum, const Matrix& mat);
        
        /// @brief Access elements
        ///
        /// Returns a pointer  to  UniTensor elements.
        /// @return Pointer to UniTensor elements
        double* getElem();
        
        
        /// @brief Copy elements
        ///
        /// Copies the first elemNum() elements from \c elem, replacing the original ones.
        /// @param elem elements to be copied.
        void setElem(const double* elem, bool _ongpu = false);
        /// @overload
        void setElem(const std::vector<double>& elem, bool _ongpu = false);
        
        /// @brief Access individual element
        ///
        /// Returns the element at linear position \c idx. The first bond’s dimension is the most significant
        ////dimension.
        ///
        /// @param    idx linear position of element
        /// @return Element at index \c idx.
        double operator[](size_t idx);
        
        /// @brief Assign elements
        ///
        /// Set all  elements to zero.
        void set_zero();
        
        /// @brief Assign elements
        ///
        /// Set all  elements in the block with quantum number \c qnum to zero.
        /// @param qnum Block quantum number
        void set_zero(const Qnum& qnum);
        
        /// @brief Assign elements
        ///
        /// Set diagonal elements of blocks to one.
        void identity();
        /// @brief Assign elements
        ///
        /// Set diagonal elements of block with \c qnum to one.
        /// @param qnum Block quantum number
        void identity(const Qnum& qnum);
        
        /// @brief Assign elements
        ///
        /// Assigns random numbers in [0, 1) to the elements.
        void randomize();
        
        /// @brief Assign elements
        ///
        /// Assigns randomly generated orthogonal bases to the elements of UniTensor.
        ///
        ///  \c Nr = row() and \c Nc = col().
        ///
        /// If the <tt> Nr < Nc </tt>, randomly generates \c Nr orthogonal basis row vectors of dimension \c Nc.
        /// If the <tt> Nr > Nc </tt>, randomly generates \c Nc orthogonal basis column vectors of
        /// dimension \c Nr.
        void orthoRand();
        
        /// @brief Assign elements
        ///
        /// Assigns randomly generated orthogonal bases to the elements of \c qnum block.
        /// @param qnum Block quantum number
        void orthoRand(const Qnum& qnum);
        
        /// @brief Access name
        ///
        /// Return the name of the UniTensor.
        std::string getName();
        
        /// @brief Assign name
        ///
        /// Assigns name to the UniTensor.
        /// @param name Name to be assigned
        void setName(const std::string& _name);
        
        /// @brief Save UniTensor to file
        ///
        /// Saves UniTensor to a file named \c fname.
        /// @param fname filename
        void save(const std::string& fname);
        
        /// @brief Permute the order of bonds
        ///
        /// Permutes the order of bonds to the order according to \c newLabels with \c inBondNum incoming bonds.
        /// @param newLabels list of new labels
        /// @param inBondNum Number of incoming bonds after permutation
        UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
        /// @overload
        UniTensor& permute(int* newLabels, int inBondNum);
        /// @brief Permute the order of bonds
        ///
        /// Rearranges the number of incoming and outgoing bonds without changing the order of the bonds.
        /// It assigns the first \c inBondNum bonds as incoming bonds and leaving the remaining bonds as
        /// outgoing bonds
        /// @param inBondNum Number of incoming bonds after permutation
        UniTensor& permute(int inBondNum);
        
        /// @brief Transpose  block elements
        ///
        /// Transpose each quantum number block. The bonds are changed from incoming to outcoming and vice versa
        /// without changing the quantum numbers on the bonds.
        UniTensor& transpose();
        
        /// @brief Combine bonds
        ///
        /// Combines  bonds with labels in \c combined_labels.
        /// The resulting bond has the same label and bondType as the bond with the first label.
        /// @param combined_labels labels to be combined
        UniTensor& combineBond(const std::vector<int>& combined_labels);
        
        /// @brief Partial trace
        ///
        /// Traces out  two bond with labels \c la and \c lb.
        /// @param la,lb Labels of the bonds to be traced out
        UniTensor& partialTrace(int la, int lb);
        
        /// @brief Trace
        /// Traces all bonds and returns the trace value.
        /// @return Trace of UniTensor
        double trace()const;
        std::vector<_Swap> exSwap(const UniTensor& Tb)const;
        void addGate(const std::vector<_Swap>& swaps);
        
        /// @brief Multiply UniTensor by a scalar and assign
        ///
        /// Performs element-wise multiplication with a scalar \c a.
        /// @param a A scalar
        UniTensor& operator*= (double a);
        
        /// @brief   Contract UniTensor with a second tensor and assign
        ///
        /// Performs tensor contraction with another UniTensor \c Tb. It contracts out the bonds of the same
        /// labels in the UniTensor and \c Tb
        /// @param Tb A second UniTensor
        UniTensor& operator*= (const UniTensor& Tb);
        
        /// @brief   Perform  element-wise addition and assign
        ///
        /// Performs element-wise addition. The tensor \c Tb to be added must be \ref{similar} to  UniTensor.
        /// @see UniTensor::similar()
        /// @param Tb A second UniTensor
        UniTensor& operator+= (const UniTensor& Tb);
        
        /// @brief Test if two tensors are similar
        ///
        /// Two UniTensor's are  similar if the bonds of two tensors are the same.
        /// @param Tb A second UniTensor
        /// @return \c True if  UniTensor is similar to \c Tb, \c false otherwise.
        bool similar(const UniTensor& Tb)const;
        
        /// @brief Test if the elements of two tensors are the same
        ///
        /// Performs the element-wise comparison of UniTensor and \c Tb. Returns \c true if all the elements are
        ///the same.
        /// @param uT  UniTensor to be compared
        /// @return \c True if the elements of  UniTensor is the same as in \c Tb, \c false otherwise.
        bool elemCmp(const UniTensor& UniT)const;
        
        /// @brief Print out raw elements
        /// Prints out raw elements of UniTensor as the following example,
        /// @code
        ///     2,0    1,0    0,0    1,0    0,0   -1,0    0,0   -1,0   -2,0
        ///    -----------------------------------------------------------------
        ///    |
        /// 2,0|  0.142  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
        ///    |
        /// 1,0|  0.000  0.952  0.000  0.916  0.000  0.000  0.000  0.000  0.000
        ///    |
        /// 0,0|  0.000  0.000  0.198  0.000  0.335  0.000  0.768  0.000  0.000
        ///    |
        /// 1,0|  0.000  0.636  0.000  0.717  0.000  0.000  0.000  0.000  0.000
        ///    |
        /// 0,0|  0.000  0.000  0.278  0.000  0.554  0.000  0.477  0.000  0.000
        ///    |
        ///-1,0|  0.000  0.000  0.000  0.000  0.000  0.394  0.000  0.783  0.000
        ///    |
        /// 0,0|  0.000  0.000  0.629  0.000  0.365  0.000  0.513  0.000  0.000
        ///    |
        ///-1,0|  0.000  0.000  0.000  0.000  0.000  0.798  0.000  0.912  0.000
        ///    |
        ///-2,0|  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.840
        ///@endcode
        ///
        ///In the above example,  UniTensor has two incoming  and two outgoing bonds, with each bond having
        /// states with Qnum's <tt>[q(1), q(0), q(-1)]</tt>. The raw elements form a 9 by 9 Matrix.
        /// The first row shows the Qnum's, \c U1 and \c parity, in the columns below and the first column
        /// shows the quantum numbers of the rows on the right.
        std::string printRawElem(bool print=true)const;
        
        /// @brief Print out  memory usage of the existing UniTensor's
        /// 
        /// Prints out the memory usage as (for example):
        ///
        /// @code{.mat}
        /// Existing Tensors: 30
        /// Allocated Elem: 2240
        /// Max Allocated Elem: 4295
        /// Max Allocated Elem for a Tensor: 924
        /// @endcode
        ///
        /// In the above example, currently there are 30 tensors and total number of existing elements is 2240.
        /// The maximum element number for now is 4295 and the maximum element number of a tensor is 924.

        static std::string profile(bool print = true);
        
        
        /// @brief Perform contraction of UniTensor
        ///
        /// Performs tensor contraction of \c Ta and \c Tb. It contracts out the bonds of the same labels
        /// in \c Ta and \c Tb.
        
        /// @note  In contrast to \ref{ operator*} as in <tt>Ta * Tb </tt>,  this function performs
        /// contraction without copying \c Ta and \c Tb. Thus it uses less memory. When the flag \c fast is set
        /// \c true, the two tensors \c Ta and \c Tb are contracted without permuting back to origin labels.
        /// @return Ta,Tb Tensors to be contracted.
        /// @param fast A flag to set if permuted back to origin labels.  If \c true, two tensor are not
        /// permuted back. Defaults to \c false
        friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
        
        /// @brief Tensor product of two tensors
        ///
        /// Performs tensor product of \c Ta and \c Tb.
        /// @param Ta,Tb Tensors to perform tensor product.
        friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
        
        /// @brief Tensor contraction
        ///
        /// Performs tensor contraction, <tt> Ta * Tb </t>. It contracts  the bonds of the same labels in \c Ta
        /// and \c Tb by making copies and then call function \ref{contract}.
        friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
        /// @brief Muliplication
        ///
        /// Performs element-wise multiplication with a scalar \c a
        /// @param a Scalar to be multiplied to UniTensor
        friend UniTensor operator* (const UniTensor& Ta, double a);
        /// @overload
        friend UniTensor operator* (double a, const UniTensor& Ta) {
            return Ta * a;
        };
        
        /// @brief Perform  element-wise addition
        /// Performs element-wise addition.
        /// @param Ta,Tb Tensors to be added
        friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
        
        /// @brief Print out UniTensor
        ///
        /// Prints out a UniTensor \c uT as(for example):
        /// @code
        ///**************** Demo ****************
        ///     ____________
        ///    |            |
        ///0___|3          3|___2
        ///    |            |
        ///1___|3          3|___3
        ///    |            |
        ///    |____________|
        ///
        ///================BONDS===============
        ///IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///IN : (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///OUT: (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///OUT: (U1 = 1, P = 0, 0)|1, (U1 = 0, P = 0, 0)|1, (U1 = -1, P = 0, 0)|1, Dim = 3
        ///
        ///===============BLOCKS===============
        ///--- (U1 = -2, P = 0, 0): 1 x 1 = 1
        ///
        ///0.840
        ///
        ///--- (U1 = -1, P = 0, 0): 2 x 2 = 4
        ///
        ///0.394  0.783
        ///
        ///0.798  0.912
        ///
        ///--- (U1 = 0, P = 0, 0): 3 x 3 = 9
        ///
        ///0.198  0.335  0.768
        ///
        ///0.278  0.554  0.477
        ///
        ///0.629  0.365  0.513
        ///
        ///--- (U1 = 1, P = 0, 0): 2 x 2 = 4
        ///
        ///0.952  0.916
        ///
        ///0.636  0.717
        ///
        ///--- (U1 = 2, P = 0, 0): 1 x 1 = 1
        ///
        ///0.142
        ///
        ///Total elemNum: 19
        ///***************** END ****************
        /// @endcode
        ///  In the above example, \c uT has four bonds with default labels [0, 1, 2, 3]. The bonds 0 and 1 are
        /// incoming bonds, and  2, 3 are out-going bonds. Each bond has three states
        /// corresponding to three U1 quantum number [-1, 0, 1]. The block elements of the
        /// tensor are als shown. There are five blocks of various <tt>U1= [-2, -1, 0, 1, 2]</tt> and
        /// various sizes. The total element number is 19.
        friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);
        
        friend class Node;
        friend class Network;
        
    private:
        std::string name;
        DOUBLE *elem;       //Array of elements
        int status; //Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements
        std::vector<Bond> bonds;
        std::map<Qnum, Block> blocks;
        std::vector<int>labels;
        void packMeta();
        int RBondNum;   //Row bond number
        int RQdim;
        int CQdim;
        size_t m_elemNum;
        std::map<int, Block*> RQidx2Blk;    //Qidx to the Block
        std::map<int, size_t> QidxEnc;
        std::map<int, size_t> RQidx2Off;    //the row offset starts from the block origin of a qnum
        std::map<int, size_t> CQidx2Off;    //the col offset starts from the block origin of a qnum
        std::map<int, size_t> RQidx2Dim;
        std::map<int, size_t> CQidx2Dim;
        bool ongpu;
        static int COUNTER;
        static int64_t ELEMNUM;
        static size_t MAXELEMNUM;
        static size_t MAXELEMTEN;   //Max number of element of a tensor
        //Private Functions
        void grouping();
        void initUniT();
        static const int HAVEBOND = 1;        /**< A flag for initialization */
        static const int HAVEELEM = 2;        /**< A flag for having element assigned */
    };
    UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast = false);
    UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
};  /* namespace uni10 */
#endif /* UNITENSOR_H */
