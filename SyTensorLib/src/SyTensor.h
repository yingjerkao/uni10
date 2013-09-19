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

typedef struct{
	int b1; 
	int b2; 
}_Swap;

/**
 * @brief Class of the symmetry tensor.
 */
class SyTensor_t{
	public:
    /**
     * @brief
     * How frequent it is used: *
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
     * @warning <tt>assert( _bonds.size() > 0 )</tt>
     */
		SyTensor_t(vector<Bond_t>& _bonds, const string& _name = "");

    /**
     * @brief To construct a tensor from a given bond array and a given label array.\n
     * How frequent it is used: * *
     * @param _bonds An STL vector of object @c Bond_t.
     * @param _labels An STL interger vector, describing the labels of bonds.
     * @param _name The given name of a tensor, STL string.
     * @see File demo/SyTensor_basic.cpp
     * @note The number of bonds must be larger than one, that is, the library does not support rank 0 tensor.
     * @note Each label is 1-1 corresponding to each bond in the order of array.
     * @warning <tt>assert( _bonds.size() > 0 )</tt>
     * @warning <tt>assert( _bonds.size() == _labels.size() )</tt>
     */
		SyTensor_t(vector<Bond_t>& _bonds, vector<int>& labels, const string& _name = "");

    /**
     * @brief To construct a tensor from a given bond array and a given label array.\n
     * How frequent it is used: * *
     * @param _bonds An STL vector of object @c Bond_t.
     * @param _labels An integer array, describing the labels of bonds.
     * @param _name The given name of a tensor, STL string.
     * @see File demo/SyTensor_basic.cpp
     * @note The number of bonds must be larger than one, that is, the library does not support rank 0 tensor.
     * @note Each label is 1-1 corresponding to each bond in the order of array.
     * @warning <tt>assert( _bonds.size() > 0 )</tt>
     * @warning <tt>assert( _bonds.size() == _labels.size() )</tt>
     */
		SyTensor_t(vector<Bond_t>& _bonds, int* labels, const string& _name = "");

    /**
     * @brief A deep copy constructor.\n
     * How frequent it is used: * * *
     * @see File demo/SyTensor_basic.cpp
     */
		SyTensor_t(const SyTensor_t& SyT);

    /**
     * @brief A deep copy assignment.\n
     * How frequent it is used: * *
     * @see File demo/SyTensor_basic.cpp
     */
		SyTensor_t& operator=(const SyTensor_t& SyT);

		~SyTensor_t();

    /**
     * @brief Add labels to the Tensor.\n
     * How frequent it is used: * * *
     * @param newLabels An STL interger vector, describing the labels of bonds.
     * @see File demo/SyTensor_basic.cpp
     * @note Each added label is 1-1 corresponding to each bond in the order of array.
     * @warning <tt>assert( _bonds.size() == _labels.size() )</tt>
     */
		void addLabel(vector<int>& newLabels);

    /**
     * @brief Add labels to the Tensor.\n
     * How frequent it is used: * * *
     * @param newLabels An interger array, describing the labels of bonds.
     * @see File demo/SyTensor_basic.cpp
     * @note Each added label is 1-1 corresponding to each bond in the order of array.
     * @warning <tt>assert( _bonds.size() == _labels.size() )</tt>
     */
		void addLabel(int* newLabels);

    /**
     * @brief Add non-blocked elements to the tensor.\n
     * How frequent it is used: * * *
     * @param rawElem An array of element type of size equal to @p elemNum.
     * @see File demo/SyTensor_basic.cpp
     * @note The alignment of the given tensor elements should follow the order of the bonds.
     */
		void addRawElem(double* rawElem);

    /**
     * @brief Get the value of correspoinding array of indices.\n
     * How frequent it is used: *
     * @para idxs An STL vector of interger array, describing the indices.
     * @see File demo/SyTensor_basic.cpp
     */
		double at(vector<int>idxs)const;

    /**
     * @brief Get an array of quantum numbers of the blocks.\n
     * How frequent it is used: * *
     * @return An STL vector of type @c Qnum_t.
     * @see File demo/SyTensor_basic.cpp
     */
    vector<Qnum_t> qnums();

    /**
     * @brief Write the tensor to an output file of filename @p fname.\n
     * How frequent it is used: *
     * @para fname A STL string, describing the filename of the output file.
     * @see File demo/SyTensor_basic.cpp
     */
		void save(const string& fname);

    /**
     * @brief Reshape the element of the tensor, that is, change the order of bonds to the order of @p newLabels and also change the element alignment to the corresponding order.\n
     * How frequent it is used: * * *
     * @param newLabels An STL interger vector, describing the labels of bonds after reshape.
     * @param rowBondNum An interger, describing the number of row bonds .
     * @see File demo/SyTensor_tool.cpp
     * @note Reshape may cause change of the order of bonds, quantum numbers of blocks of the tensor and the alignment of tensor elements.
     * @warning The only difference between @p newLabels and the original @p labels is the order of the array elements.
     */
		void reshape(vector<int>& newLabels, int rowBondNum);

    /**
     * @brief Reshape the element of the tensor, that is, change the order of bonds to the order of @p newLabels and also change the element alignment to the corresponding order.\n
     * How frequent it is used: * * *
     * @param newLabels An interger array, describing the labels of bonds after reshape.
     * @param rowBondNum An interger, describing the number of row bonds .
     * @see File demo/SyTensor_tool.cpp
     * @note Reshape may cause change of the order of bonds, quantum numbers of blocks of the tensor and the alignment of tensor elements.
     * @warning The only difference between @p newLabels and the original @p labels is the order of the array elements.
     */
		void reshape(int* newLabels, int rowBondNum);

    /**
     * @brief Transpose the tensor.\n
     * How frequent it is used: * * *
     * @see File demo/SyTensor_tool.cpp
     */
		void transpose();

    /**
     * @brief Randomly give a value(0 ~ 1.0) to each element.\n
     * How frequent it is used: *
     * @see File demo/SyTensor_tool.cpp
     */
		void randomize();

		void setName(const string& _name);
		int64_t getElemNum()const{return elemNum;};

		void check();
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
		friend void printRawElem(const SyTensor_t& SyT, const string& fname="");
		friend class Node_t;
		friend class Network_t;
		void orthoRand();
		void orthoRand(const Qnum_t& qnum);
		void eye();
		void eye(const Qnum_t& qnum);
		void bzero(const Qnum_t& qnum);
		void bzero();
		vector<bool> addSwap(vector<_Swap>swaps);
		void addGate(vector<_Swap>signs);
	private:
		string name;
		int status;	//Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements
		vector<Bond_t> bonds;
		map<Qnum_t, Block_t> blocks;
		vector<int>labels;
		DOUBLE *elem;		//Array of elements
		int RBondNum;	//Row bond number
		int64_t elemNum;
		vector<Block_t*> RQidx2Blk;	//Qidx to the Block
		vector<bool> Qidx;
		vector<int> RQidx2Off;	//the row offset starts from the block origin of a qnum
		vector<int> CQidx2Off;	//the col offset starts from the block origin of a qnum
		static int COUNTER;
		static int64_t ELEMNUM;
		static int64_t MAXELEMNUM;
		static int64_t MAXELEMTEN;	//Max number of element of a tensor
		//Private Functions
		void grouping();
		void initSyT();
};
