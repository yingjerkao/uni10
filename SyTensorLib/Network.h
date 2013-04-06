#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <set>
using namespace std;
//Bond property

class Node_t{
	private:
		SyTensor_t* SyT;
		set<int> labels;	
		map<int, Bond_t> bonds;
		int elemNum;
};

class lBond_t{	//light bond
	vector<Qnum_t>Qnums;	//Quantum numbers
	vector<int>Qdegs;	//Degeneracy in each quantum sector
};

class Network_t {
	private:
		Node_t root;
};
