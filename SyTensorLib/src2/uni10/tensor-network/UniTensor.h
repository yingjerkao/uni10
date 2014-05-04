/**
 * @file UniTensor.h
 * @author Yun-Da Hsieh
 * @date 28 Aug 2013
 * @brief This is the header file for the class of symmetry tensor "UniTensor".
 *
 * @see http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * @see http://www.stack.nl/~dimitri/doxygen/commands.html
 */
#ifndef SYTENSOR_H
#define SYTENSOR_H
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
#define DOUBLE double

#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
#include <uni10/data-structure/Block.h>
#include <uni10/tensor-network/Matrix.h>

/**
 * @brief Class of the symmetry tensor.
 */
namespace uni10{
class UniTensor{
	public:
		UniTensor(double val = 1.0);
		UniTensor(const std::string& fname);
		UniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "");
		UniTensor(const std::vector<Bond>& _bonds, std::vector<int>& labels, const std::string& _name = "");
		UniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "");
		UniTensor(const UniTensor& UniT);
		UniTensor& operator=(const UniTensor& UniT);
		UniTensor& assign(const std::vector<Bond>& _bond);
		~UniTensor();
		void addLabel(const std::vector<int>& newLabels);
		void addLabel(int* newLabels);
		void addRawElem(double* rawElem);
		double at(std::vector<int>idxs)const;
		double& operator[](size_t idx);
    	std::vector<Qnum> blockQnum()const;
    	Qnum blockQnum(int idx)const;
		size_t blockNum()const;
		void save(const std::string& fname);
		std::vector<int> label()const;
		int label(int idx)const;
		std::vector<Bond> bond()const;
		Bond bond(int idx)const;
		void setName(const std::string& _name);
		std::string getName();
		size_t elemNum()const;
		size_t bondNum()const;
		int inBondNum()const;
		static void check();
		UniTensor& permute(const std::vector<int>& newLabels, int inBondNum);
		UniTensor& permute(int* newLabels, int inBondNum);
		UniTensor& transpose();
		void randomize();
		friend std::ostream& operator<< (std::ostream& os, const UniTensor& UniT);
		friend UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast);
		friend UniTensor outer(const UniTensor& Ta, const UniTensor& Tb);
		friend UniTensor operator*(const UniTensor& Ta, const UniTensor& Tb);
		UniTensor& operator*= (const UniTensor& Tb);
		friend UniTensor operator+ (const UniTensor& Ta, const UniTensor& Tb);
		UniTensor& operator+= (const UniTensor& Tb);
		friend UniTensor operator* (const UniTensor& Ta, double a);
		friend UniTensor operator* (double a, const UniTensor& Ta){return Ta * a;};
		UniTensor& operator*= (double a);
		Matrix getBlock(const Qnum& qnum, bool diag = false)const;
		void putBlock(const Qnum& qnum, const Matrix& mat);
		std::map<Qnum, Matrix> getBlocks()const;
		Matrix rawElem()const;
		void printRawElem()const;
		friend class Node;
		friend class Network;
		void orthoRand();
		void orthoRand(const Qnum& qnum);
		void eye();
		void eye(const Qnum& qnum);
		void set_zero(const Qnum& qnum);
		void set_zero();
		std::vector<_Swap> exSwap(const UniTensor& Tb)const;
		bool similar(const UniTensor& Tb)const;
		void addGate(std::vector<_Swap> swaps);
		bool elemCmp(const UniTensor& UniT)const;
		double trace()const;
		UniTensor& combineBond(const std::vector<int>& combined_labels);
		UniTensor& partialTrace(int la, int lb);
	private:
		std::string name;
		DOUBLE *elem;		//Array of elements
		int status;	//Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements
		std::vector<Bond> bonds;
		std::map<Qnum, Block> blocks;
		std::vector<int>labels;
		void packMeta();
		int RBondNum;	//Row bond number
		int RQdim;
		int CQdim;
		int64_t m_elemNum;
		std::map<int, Block*> RQidx2Blk;	//Qidx to the Block
		std::map<int, int> QidxEnc;
		std::map<int, int> RQidx2Off;	//the row offset starts from the block origin of a qnum
		std::map<int, int> CQidx2Off;	//the col offset starts from the block origin of a qnum
		std::map<int, int> RQidx2Dim;
		std::map<int, int> CQidx2Dim;
		static int COUNTER;
		static int64_t ELEMNUM;
		static int64_t MAXELEMNUM;
		static int64_t MAXELEMTEN;	//Max number of element of a tensor
		//Private Functions
		void grouping();
		void initUniT();
		static const int HAVEBOND = 1;		  /**< A flag for initialization */
		static const int HAVEELEM = 2;		  /**< A flag for having element assigned */
		Matrix printRaw(bool flag)const;
};
UniTensor contract(UniTensor& Ta, UniTensor& Tb, bool fast = false);
UniTensor outer(const UniTensor& Ta, const UniTensor& Tb);
};	/* namespace uni10 */	
#endif /* SYTENSOR_H */
