#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
using namespace std;
//Bond property
class SyTensor_t;
class Bond_t;
class Qnum_t;

class Node_t{
	public:
		Node_t();
		Node_t(SyTensor_t* Tp);
		float assess(Node_t& nd);
		friend ostream& operator<< (ostream& os, const Node_t& nd);
	private:
		SyTensor_t* T;
		vector<int> labels;	
		vector<Bond_t> bonds;
		int64_t elemNum;
		Node_t* parent;
		Node_t* left;
		Node_t* right;
};

class lBond_t{	//light bond
	vector<Qnum_t>Qnums;	//Quantum numbers
	vector<int>Qdegs;	//Degeneracy in each quantum sector
};

class Network_t {
	public:
		Network_t();
		Network_t(vector<SyTensor_t*>& tens);
		Node_t* add(SyTensor_t*);
		SyTensor_t contract();
		void optimize(int num=1);
	private:
		vector<Node_t*> tensors;
		vector<int> order;
		Node_t* root;
		Node_t* settle(Node_t*); 
};
