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
const int INIT = 1;		//initialized
const int HAVELABEL = 2;		//tensor with label
const int HAVEELEM = 4;		//tensor with elements
const int DISPOSABLE = 8;		//The elements of tensor is disposable through 'clone', 'reshapeClone', 'reshapeElem'. The elements of disposable Tensor is read-only.
const int ELEMFREED = 16;		//the memory space of elements is freed
#define DOUBLE	double
#include "Block.h"
#include "Bond.h"
#include "myLapack.h"
class Qnum_t;
class Block_t;
class Bond_t;

class SyTensor_t{
	public:
		SyTensor_t(const string& _name = "Tensor");
		SyTensor_t(vector<Bond_t>& _bonds, const string& _name = "Tensor");
		SyTensor_t(vector<Bond_t>& _bonds, vector<int>& labels, const string& _name = "Tensor");
		SyTensor_t(vector<Bond_t>& _bonds, int* labels, const string& _name = "Tensor");
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
		void setName(string _name);
		double at(vector<int>idxs)const;
		void check();
		friend ostream& operator<< (ostream& os, SyTensor_t& SyT);
		friend SyTensor_t operator* (SyTensor_t& Ta, SyTensor_t& Tb);
		void operator*= (SyTensor_t& Tb);
		friend SyTensor_t operator+ (const SyTensor_t& Ta, const SyTensor_t& Tb);
		void operator+= (const SyTensor_t& Tb);
		friend SyTensor_t operator* (const SyTensor_t& Ta, double a);
		friend SyTensor_t operator* (double a, const SyTensor_t& Ta){return Ta * a;};
		void operator*= (double a);
		Block_t& getBlock(Qnum_t qnum);
		friend void printRawElem(const SyTensor_t& SyT);
		friend class Node_t;
		friend class Network_t;
		void orthoRand();
		void orthoRand(Qnum_t qnum);
		void eye();
		void eye(Qnum_t qnum);
		void elemset(Qnum_t qnum, double* elem, int64_t num);
		void bzero(Qnum_t qnum);
		void bzero();
	private:
		string name;
		int status;			//Check initialization, 1 initialized, 3 initialized with label, 5 initialized with elements
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

