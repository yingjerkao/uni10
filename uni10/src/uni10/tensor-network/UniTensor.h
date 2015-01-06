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
    

        /// @brief Copy UniTensor content
        ///
        /// Assigns new content to the UniTensor from \c UniT, replacing the original contents
        /// @param UniT tensor to be copied
        ///
        UniTensor& operator=(const UniTensor& UniT);
        
        /// @brief Assign bonds to UniTensor
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
        
        /// @brief Access the label of Bond \c idx
        ///
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
        
        double* getElem();
        void setElem(const double* elem, bool _ongpu = false);
        void setElem(const std::vector<double>& elem, bool _ongpu = false);
        double operator[](size_t idx);
        
        void set_zero();
        void set_zero(const Qnum& qnum);
        void identity();
        void identity(const Qnum& qnum);
        void randomize();
        void orthoRand();
        void orthoRand(const Qnum& qnum);
        std::string getName();
        void setName(const std::string& _name);
        void save(const std::string& fname);
        UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
        UniTensor& permute(int* newLabels, int inBondNum);
        UniTensor& permute(int inBondNum);
        UniTensor& transpose();
        
        UniTensor& combineBond(const std::vector<int>& combined_labels);
        UniTensor& partialTrace(int la, int lb);
        double trace()const;
        std::vector<_Swap> exSwap(const UniTensor& Tb)const;
        void addGate(const std::vector<_Swap>& swaps);
        UniTensor& operator*= (double a);
        UniTensor& operator*= (const UniTensor& Tb);
        UniTensor& operator+= (const UniTensor& Tb);
        bool similar(const UniTensor& Tb)const;
        bool elemCmp(const UniTensor& UniT)const;
        std::string printRawElem(bool print=true)const;
        static std::string profile(bool print = true);
        
        friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
        friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
        friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
        friend UniTensor operator* (const UniTensor& Ta, double a);
        friend UniTensor operator* (double a, const UniTensor& Ta) {
            return Ta * a;
        };
        friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
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
