#include "SyTensor.h"
#include "Network.h"
#include <boost/algorithm/string.hpp>
Node_t::Node_t(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL), point(0){
}

Node_t::Node_t(SyTensor_t* Tp): T(Tp), elemNum(Tp->elemNum), labels(Tp->labels), bonds(Tp->bonds), name(Tp->name), parent(NULL), left(NULL), right(NULL), point(0){	
	assert(Tp->status & Tp->INIT);
	assert(Tp->status & Tp->HAVELABEL);
}

Node_t::Node_t(const Node_t& nd): T(nd.T), elemNum(nd.elemNum), labels(nd.labels), bonds(nd.bonds), parent(nd.parent), left(nd.left), right(nd.right), point(nd.point){	
}

Node_t::Node_t(std::vector<Bond_t>& _bonds, std::vector<int>& _labels): T(NULL), labels(_labels), bonds(_bonds), parent(NULL), left(NULL), right(NULL), point(0){	
	elemNum = cal_elemNum(bonds);
}

Node_t::~Node_t(){
}

void Node_t::delink(){
	parent = NULL;
	left = NULL;
	right = NULL;
	point = 0;
}

Node_t Node_t::contract(Node_t* nd){
	int AbondNum = bonds.size();
	int BbondNum = nd->bonds.size();
	std::vector<Bond_t> cBonds;
	std::vector<int> markB(BbondNum, 0);
	std::vector<int> newLabelC;
	int conBondNum = 0;
	bool match;
	for(int a = 0; a < AbondNum; a++){
		match = false;
		for(int b = 0; b < BbondNum; b++)
			if(labels[a] == nd->labels[b]){
				markB[b] = 1;
				match = true;
				conBondNum++;
				break;
			}
		if(!match){
			newLabelC.push_back(labels[a]);
			cBonds.push_back(bonds[a]);
		}
	}
	for(int b = 0; b < BbondNum; b++)
		if(markB[b] == 0){
			newLabelC.push_back(nd->labels[b]);
			cBonds.push_back(nd->bonds[b]);
		}
	int rBondNum = AbondNum - conBondNum;
	int cBondNum = BbondNum - conBondNum;
	for(int a = 0; a < rBondNum; a++)
		cBonds[a].change(BD_ROW);
	for(int a = 0; a < cBondNum; a++)
		cBonds[rBondNum + a].change(BD_COL);

	Node_t par(cBonds, newLabelC);
	return par;
}

float Node_t::metric(Node_t* nd){	//Bigger is better
	int AbondNum = bonds.size();
	int BbondNum = nd->bonds.size();
	std::vector<Bond_t> cBonds;
	std::vector<int> markB(BbondNum, 0);
	int conBondNum = 0;
	bool match;
	for(int a = 0; a < AbondNum; a++){
		match = false;
		for(int b = 0; b < BbondNum; b++)
			if(labels[a] == nd->labels[b]){
				markB[b] = 1;
				match = true;
				conBondNum++;
				break;
			}
		if(!match)
			cBonds.push_back(bonds[a]);
	}
	if(conBondNum == 0)
		return -1;
	for(int b = 0; b < BbondNum; b++)
		if(markB[b] == 0)
			cBonds.push_back(nd->bonds[b]);
	int rBondNum = AbondNum - conBondNum;
	int cBondNum = BbondNum - conBondNum;
	for(int a = 0; a < rBondNum; a++)
		cBonds[a].change(BD_ROW);
	for(int a = 0; a < cBondNum; a++)
		cBonds[rBondNum + a].change(BD_COL);
	int64_t newElemNum = cal_elemNum(cBonds);
	return float(elemNum + nd->elemNum) / newElemNum;
}

int64_t Node_t::cal_elemNum(std::vector<Bond_t>& _bonds){
	int rBondNum = 0;
	int cBondNum = 0;
	for(int b = 0; b < _bonds.size(); b++)
		if(_bonds[b].type == BD_ROW)
			rBondNum++;
		else if(_bonds[b].type == BD_COL)
			cBondNum++;
	Qnum_t qnum(0, 0);
	int dim;
	std::map<Qnum_t,int> row_QnumMdim;
	std::vector<int> row_offs(rBondNum, 0);
	if(rBondNum){
		while(1){
			qnum.set();
			dim = 1;
			for(int b = 0; b < rBondNum; b++){
				qnum = qnum * _bonds[b].Qnums[row_offs[b]];
				dim *= _bonds[b].Qdegs[row_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end())
				row_QnumMdim[qnum] += dim;
			else
				row_QnumMdim[qnum] = dim;
			int bidx;
			for(bidx = rBondNum - 1; bidx >= 0; bidx--){
				row_offs[bidx]++;
				if(row_offs[bidx] < _bonds[bidx].Qnums.size())
					break;
				else
					row_offs[bidx] = 0;
			}
			if(bidx < 0)	//run over all row_bond offsets
				break;
		}
	}
	else{
		qnum.set();
		row_QnumMdim[qnum] = 1;
	}

	std::map<Qnum_t,int> col_QnumMdim;
	std::vector<int> col_offs(cBondNum, 0);
	if(cBondNum){
		while(1){
			qnum.set();
			dim = 1;
			for(int b = 0; b < cBondNum; b++){
				qnum = qnum * _bonds[b + rBondNum].Qnums[col_offs[b]];
				dim *= _bonds[b + rBondNum].Qdegs[col_offs[b]];
			}
			if(row_QnumMdim.find(qnum) != row_QnumMdim.end()){
				if(col_QnumMdim.find(qnum) != col_QnumMdim.end())
					col_QnumMdim[qnum] += dim;
				else
					col_QnumMdim[qnum] = dim;
			}
			int bidx;
			for(bidx = cBondNum - 1; bidx >= 0; bidx--){
				col_offs[bidx]++;
				if(col_offs[bidx] < _bonds[bidx + rBondNum].Qnums.size())
					break;
				else
					col_offs[bidx] = 0;
			}
			if(bidx < 0)	//run over all row_bond offsets
				break;
		}
	}
	else{
		qnum.set();
		if(row_QnumMdim.find(qnum) != row_QnumMdim.end())
			col_QnumMdim[qnum] = 1;
	}
	int64_t _elemNum = 0;
	std::map<Qnum_t,int>::iterator it;
	std::map<Qnum_t,int>::iterator it2;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		_elemNum += it->second * it2->second;
	}
	return _elemNum;
}

/*
Network_t::Network_t(): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
}
*/
/*
Network_t::Network_t(std::vector<SyTensor_t*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	for(int i = 0; i < tens.size(); i++)
		add(tens[i]);
}
*/

Network_t::Network_t(const std::string& fname): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	fromfile(fname);
	int Tnum = label_arr.size() - 1;
	swapflags.assign(Tnum, false);
	std::vector<_Swap> swaps;
	swaps_arr.assign(Tnum, swaps);
	leafs.assign(Tnum, NULL);
	tensors.assign(Tnum, NULL);
}

Network_t::Network_t(const std::string& fname, const std::vector<SyTensor_t*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	fromfile(fname);
	assert((label_arr.size() - 1) == tens.size());
	int Tnum = tens.size();
	swapflags.assign(Tnum, false);
	std::vector<_Swap> swaps;
	swaps_arr.assign(Tnum, swaps);
	leafs.assign(Tnum, NULL);
	tensors.assign(Tnum, NULL);
	for(int i = 0; i < Tnum; i++){
		if(tens[i]->name.length() > 0)
			assert(tens[i]->name == names[i]);
		assert(tens[i]->RBondNum == Rnums[i]);
		SyTensor_t* ten = new SyTensor_t(*(tens[i]));
		ten->setName(names[i]);
		ten->addLabel(label_arr[i]);
		tensors[i] = ten;
		Node_t* ndp = new Node_t(ten);
		leafs[i] = ndp;
	}
}

void Network_t::fromfile(const std::string& fname){
	std::string str;
	std::ifstream infile;
	infile.open (fname.c_str());
	int lnum = 0;
	int MAXLINES = 1000;
	int pos = 0;
	int endpos = 0;
	std::string tar("1234567890-");
	std::vector<std::string> ord;
	while(lnum < MAXLINES){
		getline(infile, str); // Saves the line in STRING.
		if(infile.eof())
			break;
		pos = str.find(":");
		if(pos == std::string::npos)
			break;
		std::string name = str.substr(0, pos);
		boost::algorithm::trim(name);
		if(name == "ORDER"){
			std::string del(" ,;");
			while(((pos = str.find_first_not_of(del, pos + 1)) != std::string::npos)){
				endpos = str.find_first_of(del, pos + 1);
				if(endpos == std::string::npos)
					ord.push_back(str.substr(pos));
				else
					ord.push_back(str.substr(pos, endpos - pos));
				pos = endpos;
				if(pos == std::string::npos)
					break;
			}
			break;
		}
		names.push_back(name);
		std::vector<int> labels;
		int Rnum = 0;
		int cnt = 0;
		int tmp;
		while((pos = str.find_first_of(tar, pos + 1)) != std::string::npos){
			if(Rnum == 0){
				tmp = str.find(";", endpos);
				if(tmp != std::string::npos && tmp < pos)
					Rnum = cnt;
			}
			endpos = str.find_first_not_of(tar, pos + 1);
			std::string label;
			if(endpos == std::string::npos)
				label = str.substr(pos);
			else
				label = str.substr(pos, endpos - pos);
			char* pEnd;
			labels.push_back(strtol(label.c_str(), &pEnd, 10));
			pos = endpos;
			if(Rnum == 0)
				cnt++;
			if(pos == std::string::npos)
				break;
		}
		label_arr.push_back(labels);
		Rnums.push_back(Rnum);
		lnum ++;
	}
	int numT = names.size() - 1;
	assert(names[numT] == "TOUT");
	assert(names.size() > 0);
	order.assign(numT, 0);
	if(ord.size() == numT){
		bool found;
		for(int i = 0; i < numT; i++){
			found = false;
			for(int j = 0; j < numT; j++){
				if(ord[i] == names[j]){
					found = true;
					order[i] = j;
					break;
				}
			}
			assert(found);
		}
	}
	else{
		for(int i = 0; i < numT; i++)
			order[i] = i;
	}
	infile.close();
}

/*
Node_t* Network_t::add(SyTensor_t* SyT){
	assert(label_arr.size() == 0);
	Node_t* ndp = new Node_t(SyT);
	order.push_back(leafs.size());
	leafs.push_back(ndp);
	return ndp;
}
*/

Node_t* Network_t::replaceWith(int idx, SyTensor_t* SyT, bool force){
	assert(label_arr.size() > 0 && idx >= 0 && idx < (label_arr.size()-1));
	if((!force) && load)
		destruct();
	if(SyT->name.length() > 0)
		assert(SyT->name == names[idx]);
	assert(SyT->RBondNum == Rnums[idx]);

	if(leafs[idx] != NULL){
		assert(tensors[idx]->similar(*SyT));
		*(tensors[idx]) = *SyT;
		tensors[idx]->addLabel(label_arr[idx]);
		tensors[idx]->setName(names[idx]);
		swapflags[idx] = false;
	}
	else{
		SyTensor_t* ten = new SyTensor_t(*SyT);
		ten->setName(names[idx]);
		ten->addLabel(label_arr[idx]);
		tensors[idx] = ten;
		Node_t* ndp = new Node_t(ten);
		leafs[idx] = ndp;
	}
	return leafs[idx];
}

void Network_t::branch(Node_t* sbj, Node_t* tar){
	Node_t* par = new Node_t(tar->contract(sbj));
	if(sbj->parent == NULL){	//create a parent node
		if(tar->parent != NULL){	//tar is not root
			if(tar->parent->left == tar)	// tar on the left side of its parent
				tar->parent->left = par;
			else
				tar->parent->right = par;
			par->parent = tar->parent;
		}
		else{	//tar is root
			par->parent = NULL;
			root = par;
		}
	}
	else{	//sbj and tar have same parent and replace the parent node
		if(tar->parent->parent != NULL){
			if(tar->parent->parent->left == tar->parent)	// tar on the left side of its parent
				tar->parent->parent->left = par;
			else
				tar->parent->parent->right = par;
			par->parent = tar->parent->parent;
		}
		else{	//tar->parent is root
			par->parent = NULL;
			root = par;
		}
		delete tar->parent;
	}
	par->left = tar;
	par->right = sbj;
	tar->parent = par;
	sbj->parent = par;
	par->point = tar->metric(sbj);

	if(sbj->parent->parent != NULL){	//propagate up
		sbj = sbj->parent;
		branch(sbj->parent->right, sbj->parent->left);
	}
}

void Network_t::matching(Node_t* sbj, Node_t* tar){
	if(tar == NULL){	//tar is root
		root = sbj;
	}
	else if(tar->T == NULL){	//not leaf
		if(sbj->metric(tar) > 0){	//has contracted bonds
			assert(tar->left != NULL && tar->right != NULL);
			float tar_p = tar->point;
			float lft_p = 0, rht_p = 0;
			if((lft_p = sbj->metric(tar->left)) > tar_p || (rht_p = sbj->metric(tar->right)) > tar_p){	//go deeper comparison to the children
				if(lft_p > rht_p)
					matching(sbj, tar->left);
				else
					matching(sbj, tar->right);
			}
			else	//contract!!!
				branch(sbj, tar);
		}
		else	//contract!!!
			branch(sbj, tar);
	}
	else{	//contract!!!
		branch(sbj, tar);
	}
}

void Network_t::clean(Node_t* nd){
	if(nd->T != NULL)	//leaf
		return;
	clean(nd->left);
	clean(nd->right);
	delete nd;
}

void Network_t::destruct(){
	clean(root);
	root = NULL;
	for(int i = 0; i < leafs.size(); i++)
		leafs[i]->delink();
	conOrder.clear();
	for(int t = 0; t < tensors.size(); t++){
		if(swapflags[t]){
			tensors[t]->addGate(swaps_arr[t]);
			swapflags[t] = false;
		}
		swaps_arr[t].clear();
	}
	load = false;
}

void Network_t::construct(){
	for(int i = 0; i < order.size(); i++){
		assert(leafs[order[i]] != NULL);
		matching(leafs[order[i]], root);
	}
	addSwap();
	load = true;
}

//void Network_t::optimize(int num){
//	assert(false);//not a ready function
//}

SyTensor_t Network_t::launch(const std::string& _name){
	if(!load)
		construct();
	for(int t = 0; t < tensors.size(); t++)
		if(!swapflags[t]){
			tensors[t]->addGate(swaps_arr[t]);
			swapflags[t] = true;
		}
	SyTensor_t SyT = merge(root);
	int idx = label_arr.size() - 1;
	if(label_arr.size() > 0)
		SyT.reshape(label_arr[idx], Rnums[idx]);
	SyT.setName(_name);
	return SyT;
}

SyTensor_t Network_t::merge(Node_t* nd){
	if(nd->left->T == NULL){
		SyTensor_t lftT = merge(nd->left);
		if(nd->right->T == NULL){
			SyTensor_t rhtT = merge(nd->right);
			return lftT * rhtT;
		}
		else{
			return lftT * *(nd->right->T);
		}
	}
	else{
		if(nd->right->T == NULL){
			SyTensor_t rhtT = merge(nd->right);
			return *(nd->left->T) * rhtT;
		}
		else{
			return *(nd->left->T) * *(nd->right->T);
		}
	}
}

Network_t::~Network_t(){
	if(load)
		destruct();
	for(int i = 0; i < leafs.size(); i++)
		delete leafs[i];
	for(int i = 0; i < tensors.size(); i++)
		delete tensors[i];
}

void Network_t::preprint(std::ostream& os, Node_t* nd, int layer){
	if(nd == NULL)
		return;
	for(int i = 0; i < layer; i++)
		os<<"|   ";
	if(nd->T)
		os<<nd->name << "(" << nd->elemNum << "): ";
	else
		os<<"*("<<nd->elemNum<<"): ";
	for(int i = 0; i < nd->labels.size(); i++)
		os<< nd->labels[i] << ", ";
	os<<std::endl;
	preprint(os, nd->left, layer+1);
	preprint(os, nd->right, layer+1);
}

std::ostream& operator<< (std::ostream& os, Network_t& net){
	if(!net.load)
		net.construct();
	net.preprint(os, net.root, 0);
	return os;
}
std::ostream& operator<< (std::ostream& os, const Node_t& nd){
	os << "Tensor: " << nd.T<<std::endl;
	os << "elemNum: " << nd.elemNum<<std::endl;
	os << "parent: " << nd.parent<<std::endl;
	os << "left: " << nd.left<<std::endl;
	os << "right: " << nd.right<<std::endl;
	os << "labels: ";
	for(int i = 0; i < nd.labels.size(); i++)
		os << nd.labels[i] << ", ";
	os << std::endl;
	for(int i = 0; i < nd.bonds.size(); i++)
		os << "    " <<  nd.bonds[i];
	return os;
}

void Network_t::findConOrd(Node_t* nd){
	if(nd == NULL || conOrder.size() == tensors.size())
		return;
	if(nd->T){
		bool found = false;
		for(int i = 0; i < tensors.size(); i++)
			if(nd->T == tensors[i]){
				conOrder.push_back(i);
				found = true;
				break;
			}
		assert(found);
	}
	findConOrd(nd->left);
	findConOrd(nd->right);
}

void Network_t::addSwap(){
	int Tnum = leafs.size();
	findConOrd(root);
	assert(Tnum == conOrder.size());
	int tenOrder[conOrder.size()];
	memcpy(tenOrder, &(conOrder[0]), Tnum * sizeof(int));
	std::vector<_Swap> tenSwaps = _recSwap(tenOrder, Tnum);
	std::vector<_Swap> swtmp;
	for(int s = 0; s < tenSwaps.size(); s++){
		swtmp = tensors[tenSwaps[s].b1]->exSwap(*(tensors[tenSwaps[s].b2]));
		swaps_arr[tenSwaps[s].b1].insert(swaps_arr[tenSwaps[s].b1].end(), swtmp.begin(), swtmp.end());
	}
	//Distinct Swaps of each tensors
	for(int t = 0; t < Tnum; t++){
		std::map<int, bool> recs;
		std::map<int, bool>::iterator it;
		int bondNum = tensors[t]->bonds.size();
		int is, ib;
		int hash;
		for(int s = 0; s < swaps_arr[t].size(); s++){
			if(swaps_arr[t][s].b1 < swaps_arr[t][s].b2){
				is = swaps_arr[t][s].b1;
				ib = swaps_arr[t][s].b2;
			}
			else{
				is = swaps_arr[t][s].b2;
				ib = swaps_arr[t][s].b1;
			}
			int hash = is * bondNum + ib;
			if((it = recs.find(hash)) != recs.end())
				it->second ^= true;
			else
				recs[hash] = true;
		}
		swaps_arr[t].clear();
		_Swap sp;
		for (it=recs.begin(); it!=recs.end(); it++){
			if(it->second){
				sp.b1 = (it->first) / bondNum;
				sp.b2 = (it->first) % bondNum;
				swaps_arr[t].push_back(sp);
			}
		}
	}
}
