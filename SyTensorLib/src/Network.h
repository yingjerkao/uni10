#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
//Bond property
class SyTensor_t;
class Bond_t;
class Qnum_t;
class Node_t{
	public:
		Node_t();
		Node_t(SyTensor_t* Tp);
		Node_t(const Node_t& nd);
		Node_t(std::vector<Bond_t>& _bonds, std::vector<int>& _labels);
		~Node_t();
		Node_t contract(Node_t* nd);
		float metric(Node_t* nd);
		friend std::ostream& operator<< (std::ostream& os, const Node_t& nd);
		friend class Network_t;
	private:
		SyTensor_t* T;	//if T != NULL, it is leaf node
		std::vector<int> labels;
		std::vector<Bond_t> bonds;
		int64_t elemNum;
		std::string name;
		Node_t* parent;
		Node_t* left;
		Node_t* right;
		float point;
		int64_t cal_elemNum(std::vector<Bond_t>& _bonds);
		void delink();
};

class Network_t {
	public:
		//Network_t();
		Network_t(std::vector<SyTensor_t*>& tens);
		Network_t(const std::string& fname, const std::vector<SyTensor_t*>& tens);
		Network_t(std::vector< std::vector<int> > _label_arr);
		Network_t(const std::string& fname);
		~Network_t();
		//Node_t* add(SyTensor_t*);
		Node_t* replaceWith(int idx, SyTensor_t* SyT, bool force=false);	//if force is true, force replace without change the all network
		SyTensor_t launch(const std::string& name="");
		SyTensor_t launch(int* outLabels, int Rnum = 0, const std::string& name="");
		SyTensor_t launch(std::vector<int>& outLabels, int Rnum = 0, const std::string& name="");
		//void optimize(int num=1);
		friend std::ostream& operator<< (std::ostream& os, Network_t& nd);
	private:
		void preprint(std::ostream& os, Node_t* nd, int layer);	//pre-order print
		std::vector<std::string> names;
		std::vector< std::vector<int> > label_arr;
		std::vector< int > Rnums;
		std::vector<Node_t*> leafs;
		std::vector<SyTensor_t*> tensors;
		std::vector< std::vector<_Swap> > swaps_arr;
		std::vector<bool> swapflags;
		std::vector<int> conOrder;	//contraction order;
		std::vector<int> order;	//add order
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
		void fromfile(const std::string& fname);
		void findConOrd(Node_t* nd);
		void addSwap();
};
