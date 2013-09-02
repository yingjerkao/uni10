/**
 * @file SyTensor.h
 * @author Yun-Da Hsieh
 * @date 28 Aug 2013
 * @brief This is the header file for the class of symmetry tensor "SyTensor_t".
 *
 * @see http://www.stack.nl/~dimitri/doxygen/docblocks.html
 * @see http://www.stack.nl/~dimitri/doxygen/commands.html
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <assert.h>
#include <stdint.h>
using namespace std;
#include "Block.h"
#include "Bond.h"
#include "myLapack.h"
#include "Matrix.h"
#define DOUBLE	double

const int INIT = 1;		      /**< A flag for initialization */
const int HAVELABEL = 2;		/**< A flag for having labels added */
const int HAVEELEM = 4;		  /**< A flag for having element assigned */


/**
 * @brief Class of the symmetry tensor.
 */
class SyTensor_t{
	public:
    /**
     * @brief
     * How frequent it is used: *
     * @verbatim How frequent it is used: * @endverbatim
     * @see File demo/SyTensor_basic.cpp
     */
		SyTensor_t();
    /**
     * @brief To read in a binary file of a tensor which is written out by member function @c save().\n
     * How frequent it is used: * * *
     * @param fname The file name of the tensor being loaded, which is of type STL @c string.
     * @see File demo/SyTensor_basic.cpp
     */
		SyTensor_t(const string& fname);
    /**
     * @brief To construct a tensor from a given bond array.\n
     * How frequent it is used: * * *
     * @param _bonds an STL vector of object @c Bond_t.
     * @param _name The given name of a tensor, STL string.
     * @see File demo/SyTensor_basic.cpp
     * @note The number of bonds must be larger than one, that is, the library does not support rank 0 tensor.
     * @warning <tt>assert(_bonds.size() > 0)</tt>
     */
		SyTensor_t(vector<Bond_t>& _bonds, const string& _name = "");
    /**
     * @brief To construct a tensor from a given bond array and a given label array.\n
     * How frequent it is used: * *
     * @param _bonds an STL vector of object @c Bond_t.
     * @param _labels.
     * @see File demo/SyTensor_basic.cpp
     * @note each label is 1-1 corresponding to each bond in order of array.
     * @warning <tt>assert(_bonds.size() == _labels.size())</tt>
     */
		SyTensor_t(vector<Bond_t>& _bonds, vector<int>& labels, const string& _name = "");
		SyTensor_t(vector<Bond_t>& _bonds, int* labels, const string& _name = "");
		SyTensor_t(const SyTensor_t& SyT);
		~SyTensor_t();
		SyTensor_t& operator=(const SyTensor_t& SyT);
		void addLabel(vector<int>& newLabels);
		void addLabel(int* newLabels);
		void reshape(vector<int>& newLabels, int rowBondNum);
		void reshape(int* newLabels, int rowBondNum);
		void addRawElem(double* rawElem);
		void transpose();
		void randomize();
		vector<Qnum_t> qnums();
		void setName(const string& _name);
		double at(vector<int>idxs)const;
		void check();
		void save(const string& fname);
		friend ostream& operator<< (ostream& os, SyTensor_t& SyT);
		friend SyTensor_t operator* (SyTensor_t& Ta, SyTensor_t& Tb);
		void operator*= (SyTensor_t& Tb);
		friend SyTensor_t operator+ (const SyTensor_t& Ta, const SyTensor_t& Tb);
		void operator+= (const SyTensor_t& Tb);
		friend SyTensor_t operator* (const SyTensor_t& Ta, double a);
		friend SyTensor_t operator* (double a, const SyTensor_t& Ta){return Ta * a;};
		void operator*= (double a);
		Matrix_t getBlock(Qnum_t qnum, bool diag = false);
		void putBlock(const Qnum_t& qnum, Matrix_t& mat);
		friend void printRawElem(const SyTensor_t& SyT);
		friend class Node_t;
		friend class Network_t;
		void orthoRand();
		void orthoRand(const Qnum_t& qnum);
		void eye();
		void eye(const Qnum_t& qnum);
		void bzero(const Qnum_t& qnum);
		void bzero();
	private:
		string name;
		int status;	//Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements
		vector<Bond_t> bonds;
		map<Qnum_t, Block_t> blocks;
		vector<int>labels;
		DOUBLE *elem;		//Array of elements
		int RBondNum;	//Row bond number
		int64_t elemNum;
		vector<Block_t*> RQidx2Blk;
		vector<bool> Qidx;
		vector<int> RQidx2Off;
		vector<int> CQidx2Off;
		static int COUNTER;
		static int64_t ELEMNUM;
		static int64_t MAXELEMNUM;
		static int64_t MAXELEMTEN;	//Max number of element of a tensor
		//Private Functions
		void grouping();
		void initSyT();
};
