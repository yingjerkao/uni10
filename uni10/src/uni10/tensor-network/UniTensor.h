/****************************************************************************
 *  @file CMakeLists.txt
 *  @license
 *   Universal Tensor Network Library
 *   Copyright (c) 2013-2014
 *   National Taiwan University
 *   National Tsing-Hua University
 *
 *   This file is part of Uni10, the Universal Tensor Network Library.
 *
 *   Uni10 is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Uni10 is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public License
 *   along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
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
#include <cmath>
#include <cstdio>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <assert.h>
#include <sstream>
#include <stdexcept>
#include <uni10/datatype.hpp>
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/data-structure/Block.h>
#include <uni10/tensor-network/Matrix.h>
/// @brief Uni10 - the Universal Tensor %Network Library


namespace uni10 {

    class Block;
    class CBlock;
    class Matrix;
    class CMatrix;
    class CUniTensor;
    class UniTensor {

    public:
        
        /*********************  OPERATOR **************************/	    
        
        friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);
       
       
        /********************* going move **************************/	    
       
        UniTensor(muType tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(muType tp, const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(muType tp, const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
      
        /*********************  NO TYPE **************************/	    
     
        UniTensor();
        UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        UniTensor(const UniTensor& UniT);
        UniTensor(const std::string& fname);
        
        UniTensor(const Block& UniT);
        ~UniTensor();
        
        int typeID()const; 
        void setLabel(const std::vector<int>& newLabels);
        void setLabel(int* newLabels);
        std::vector<int> label()const;
        int label(size_t idx)const;
        std::string getName();
        void setName(const std::string& _name);
        size_t bondNum()const;
        size_t inBondNum()const;
        std::vector<Bond> bond()const;
        Bond bond(size_t idx)const;
        size_t elemNum()const;
        size_t blockNum()const;
        std::vector<Qnum> blockQnum()const;
        Qnum blockQnum(size_t idx)const;
        const std::map<Qnum, Block>& const_getBlocks()const;
        const Block& const_getBlock()const;
        const Block& const_getBlock(const Qnum& qnum)const;
        std::map<Qnum, Matrix> getBlocks()const;
        Matrix getBlock(bool diag = false)const;
        Matrix getBlock(const Qnum& qnum, bool diag = false)const;
        void set_zero();
        void set_zero(const Qnum& qnum);
        void identity();
        void identity(const Qnum& qnum);
        void randomize();
        void orthoRand();
        void orthoRand(const Qnum& qnum);
        void save(const std::string& fname);
        UniTensor& transpose();
        
        /*********************  REAL **********************/
        
        UniTensor(double val);
        UniTensor(rflag _tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(rflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(rflag _tp, const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        void setRawElem(const std::vector<double>& rawElem);
        void setRawElem(const double* rawElem);
        void setElem(const double* elem, bool _ongpu = false);
        void setElem(const std::vector<double>& elem, bool _ongpu = false);
        std::map<Qnum, Matrix> getBlocks(rflag _tp)const;
        Matrix getBlock(rflag _tp, bool diag = false)const;
        Matrix getBlock(rflag _tp, const Qnum& qnum, bool diag = false)const;
        void set_zero(rflag _tp);
        void set_zero(rflag _tp, const Qnum& qnum);
        void identity(rflag _tp);
        void identity(rflag _tp, const Qnum& qnum);
        void randomize(rflag _tp);
        void orthoRand(rflag _tp);
        void orthoRand(rflag _tp, const Qnum& qnum);
        UniTensor& transpose(rflag _tp);
       
        /*********************  COMPLEX **********************/
      
        UniTensor(std::complex<double> val);
        UniTensor(cflag _tp, const std::vector<Bond>& _bonds, const std::string& _name = "");
        UniTensor(cflag _tp, const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
        UniTensor(cflag _tp, const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
        void setRawElem(const std::vector< std::complex<double> >& rawElem);
        void setRawElem(const std::complex<double>* rawElem);
        void setElem(const std::complex<double>* c_elem, bool _ongpu = false);
        void setElem(const std::vector< std::complex<double> >& c_elem, bool _ongpu = false);
        std::map<Qnum, Matrix> getBlocks(cflag _tp)const;
        Matrix getBlock(cflag _tp, bool diag = false)const;
        Matrix getBlock(cflag _tp, const Qnum& qnum, bool diag = false)const;
        void set_zero(cflag _tp);
        void set_zero(cflag _tp, const Qnum& qnum);
        void identity(cflag _tp);
        void identity(cflag _tp, const Qnum& qnum);
        void randomize(cflag _tp);
        void orthoRand(cflag _tp);
        void orthoRand(cflag _tp, const Qnum& qnum);
        UniTensor& transpose(cflag _tp);
     
        /******************Friend funcs*******************/
        
        /*****************************************************/
    
        UniTensor& operator*= (double a);
        UniTensor& operator*= (std::complex<double> a);
        UniTensor& operator*= (const UniTensor& Tb);
        UniTensor& operator+= (const UniTensor& Tb);
        friend UniTensor operator*(const UniTensor& Ta, const std::complex<double>& a);
        friend UniTensor operator*(const std::complex<double>& a, const UniTensor& Ta);
        friend UniTensor operator*(const UniTensor& Ta, double a);
        friend UniTensor operator*(double a, const UniTensor& Ta);
        friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
        friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
        friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
        friend UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
        
        std::string printRawElem(bool print=true)const;
        static std::string profile(bool print = true);

        void setRawElem(const Block& blk);
        void putBlock(const Block& mat);
        void putBlock(const Qnum& qnum, const Block& mat);
        
        
        UniTensor& combineBond(const std::vector<int>& combined_labels);
        std::vector<_Swap> exSwap(const UniTensor& Tb)const;
        void addGate(const std::vector<_Swap>& swaps);
        
        bool similar(const UniTensor& Tb)const;
        bool elemCmp(const UniTensor& UniT)const;
        void clear();
        
        friend void RtoC(UniTensor& UniT);
        
        UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
        UniTensor& permute(int* newLabels, int inBondNum);
        UniTensor& permute(int inBondNum);
        
        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum = 0)const;
        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum, std::vector<Matrix>& Ls)const;
        std::vector<UniTensor> hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls)const;
        std::vector<UniTensor> hosvd(size_t modeNum, std::vector<Matrix>& Ls)const;
        std::vector<UniTensor> hosvd(size_t modeNum, std::vector<std::map<Qnum, Matrix> >& Ls)const;
       
        UniTensor& assign(muType _tp, const std::vector<Bond>& _bond);
        
        //double* getElem();
        double* getRealElem();
        std::complex<double>* getComplexElem();
        Matrix getRawElem()const;
        std::complex<double> trace()const;
        UniTensor& partialTrace(int la, int lb);
        std::complex<double> operator[](size_t idx) const;
        std::complex<double> at(const std::vector<int>& idxs)const;
        std::complex<double> at(const std::vector<size_t>& idxs)const;
        
        /****************************/

        /******** check func ********/
        bool isCelemEmpty(); 
        bool isElemEmpty();
        /****************************/
        
        /******** find bug **********/
        UniTensor& operator=(const UniTensor& UniT);
        UniTensor& assign(const std::vector<Bond>& _bond);
        /****************************/
        
        std::complex<double> at(muType _tp, const std::vector<int>& idxs)const;
        std::complex<double> at(muType _tp, const std::vector<size_t>& idxs)const;
        
        friend class CUniTensor;
        friend class Node;
        friend class Network;

    private:
        muType u_type;

        rflag r_flag;
        cflag c_flag;
        std::string name;
        double *elem;       //Array of elements
        std::complex<double>* c_elem;       //Array of elements
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
        /********************* going move **************************/	    
        void initPrototype(muType tp = EMPTY);
        void initUniT(muType tp);
        void uelemAlloc();
        void initBlocks(double* _elem);
        void initBlocks(std::complex<double>* _c_elem);
        void uelemBzero();
        /*********************  NO TYPE **************************/	    
        void initUniT();
        void initUniT(int _typeID);
        size_t grouping();
        void uelemFree();
        /*********************  REAL **********************/
        void initUniT(rflag _tp);
        size_t grouping(rflag _tp);
        void uelemAlloc(rflag _tp);
        void initBlocks(rflag _tp);
        void uelemBzero(rflag _tp);
        /*********************  COMPLEX **********************/
        void initUniT(cflag _tp);
        size_t grouping(cflag _tp);
        void uelemAlloc(cflag _tp);
        void initBlocks(cflag _tp);
        void uelemBzero(cflag _tp);
        /*****************************************************/
        
        std::vector<UniTensor> _hosvd(size_t modeNum, size_t fixedNum, std::vector<std::map<Qnum, Matrix> >& Ls, bool returnL)const;
        static const int HAVEBOND = 1;        /**< A flag for initialization */
        static const int HAVEELEM = 2;        /**< A flag for having element assigned */
    };
    UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast = false);
    UniTensor otimes(const UniTensor& Ta, const UniTensor& Tb);
    void RtoC(UniTensor& UniT);
};  /* namespace uni10 */
#endif /* UNITENSOR_H */
