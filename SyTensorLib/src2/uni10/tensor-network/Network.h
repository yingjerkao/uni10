#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
//Bond property
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
namespace uni10{
class Network {
	public:
		Network(const std::string& fname, const std::vector<UniTensor*>& tens);
		Network(const std::string& fname);
		~Network();
		//Node* add(UniTensor*);
		void putTensor(int idx, const UniTensor* UniT, bool force=false);	//if force is true, force replace without change the all network
		UniTensor launch(const std::string& name="");
		//void optimize(int num=1);
		friend std::ostream& operator<< (std::ostream& os, Network& nd);
	private:
		void preprint(std::ostream& os, Node* nd, int layer);	//pre-order print
		std::vector<std::string> names;
		std::vector< std::vector<int> > label_arr;
		std::vector< int > Rnums;
		std::vector<Node*> leafs;
		std::vector<UniTensor*> tensors;
		std::vector< std::vector<_Swap> > swaps_arr;
		std::vector<bool> swapflags;
		std::vector<int> conOrder;	//contraction order;
		std::vector<int> order;	//add order
		std::vector<int> brakets;	//add order
		Node* root;
		bool load;	//whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
		int times;	//construction times
		int tot_elem;	//total memory ussage
		int max_elem;	//maximum
		void construct();
		void destruct();
		void matching(Node* sbj, Node* tar);
		void branch(Node* sbj, Node* tar);
		UniTensor merge(Node* nd);
		void clean(Node* nd);
		void fromfile(const std::string& fname);
		void findConOrd(Node* nd);
		void addSwap();
};
};	/* namespace uni10 */	
#endif /* NETWORK_H */
