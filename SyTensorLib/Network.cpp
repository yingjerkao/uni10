#include "SyTensor.h"
#include "Network.h"
Node_t::Node_t(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL){
}
Node_t::Node_t(SyTensor_t* Tp): T(Tp), elemNum(Tp->elemNum), labels(Tp->labels), bonds(Tp->bonds), parent(NULL), left(NULL), right(NULL){	
	assert(Tp->status & INIT);
	assert(Tp->status & HAVELABEL);
}
float Node_t::assess(Node_t& nd){
	int AbondNum = bonds.size();
	int BbondNum = nd.bonds.size();
	vector<Bond_t> cBonds;
	vector<int> markB(BbondNum, 0);
	vector<int> newLabelC;
	bool match;
	for(int a = 0; a < AbondNum; a++){
		match = false;
		for(int b = 0; b < BbondNum; b++)
			if(labels[a] == nd.labels[b]){
				markB[b] = 1;
				match = true;
				break;
			}
		if(!match){
			newLabelC.push_back(labels[a]);
			cBonds.push_back(bonds[a]);
		}
	}
	for(int b = 0; b < BbondNum; b++)
		if(markB[b] == 0){
			newLabelC.push_back(nd.labels[b]);
			cBonds.push_back(nd.bonds[b]);
		}
	return 1.0f;
}
Network_t::Network_t(){
}

ostream& operator<< (ostream& os, const Node_t& nd){
	os << "Tensor: " << nd.T<<endl;
	os << "elemNum: " << nd.elemNum<<endl;
	os << "parent: " << nd.parent<<endl;
	os << "left: " << nd.left<<endl;
	os << "right: " << nd.right<<endl;
	os << "labels: ";
	for(int i = 0; i < nd.labels.size(); i++)
		os << nd.labels[i] << ", ";
	os << endl;
	for(int i = 0; i < nd.bonds.size(); i++)
		os << "    " <<  nd.bonds[i];
	return os;
}
