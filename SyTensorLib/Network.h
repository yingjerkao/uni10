#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
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
		Node_t(SyTensor_t& Tp);
		Node_t(const Node_t& nd);
		Node_t(vector<Bond_t>& _bonds, vector<int>& _labels);
		~Node_t();
		Node_t contract(Node_t* nd);
		float metric(Node_t* nd);
		friend ostream& operator<< (ostream& os, const Node_t& nd);
		friend class Network_t;
	private:
		SyTensor_t* T;
		vector<int> labels;	
		vector<Bond_t> bonds;
		int64_t elemNum;
		Node_t* parent;
		Node_t* left;
		Node_t* right;
		float point;
		int64_t cal_elemNum(vector<Bond_t>& _bonds);
};

class lBond_t{	//light bond
	vector<Qnum_t>Qnums;	//Quantum numbers
	vector<int>Qdegs;	//Degeneracy in each quantum sector
};

class Network_t {
	public:
		Network_t();
		Network_t(vector<SyTensor_t*>& tens);
		Network_t(const string& fname, vector<SyTensor_t*>& tens);
		Network_t(vector< vector<int> > _label_arr);
		Network_t(const string& fname);
		~Network_t();
		Node_t* add(SyTensor_t&);
		SyTensor_t launch(const string& name="");
		SyTensor_t launch(int* outLabels, int Rnum = 0, const string& name="");
		SyTensor_t launch(vector<int>& outLabels, int Rnum = 0, const string& name="");
		//void optimize(int num=1);
		friend ostream& operator<< (ostream& os, Network_t& nd);
	private:
		void preprint(ostream& os, Node_t* nd, int layer);	//pre-order print
		vector<string> names;
		vector< vector<int> > label_arr;
		vector<Node_t*> leafs;
		vector<int> order;
		Node_t* root;
		bool load;	//whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
		int times;	//construction times
		int tot_elem;	//total memory ussage
		int max_elem;	//maximum 
		void construct(); 
		void destruct(); 
		void matching(Node_t* sbj, Node_t* tar);
		void branch(Node_t* sbj, Node_t* tar);
		SyTensor_t merge(Node_t* nd);
		void clean(Node_t* nd);
		void fromfile(const string& fname);
};
