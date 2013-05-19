#include "SyTensor.h"
#include "Network.h"
Node_t::Node_t(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL), point(0){
}

Node_t::Node_t(SyTensor_t& Tp): T(&Tp), elemNum(Tp.elemNum), labels(Tp.labels), bonds(Tp.bonds), parent(NULL), left(NULL), right(NULL), point(0){	
	assert(Tp.status & INIT);
	assert(Tp.status & HAVELABEL);
}

Node_t::Node_t(const Node_t& nd): T(nd.T), elemNum(nd.elemNum), labels(nd.labels), bonds(nd.bonds), parent(nd.parent), left(nd.left), right(nd.right), point(nd.point){	
}

Node_t::Node_t(vector<Bond_t>& _bonds, vector<int>& _labels): T(NULL), labels(_labels), bonds(_bonds), parent(NULL), left(NULL), right(NULL), point(0){	
	elemNum = cal_elemNum(bonds);
}

Node_t::~Node_t(){
}

Node_t Node_t::contract(Node_t* nd){
	int AbondNum = bonds.size();
	int BbondNum = nd->bonds.size();
	vector<Bond_t> cBonds;
	vector<int> markB(BbondNum, 0);
	vector<int> newLabelC;
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
	vector<Bond_t> cBonds;
	vector<int> markB(BbondNum, 0);
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

int64_t Node_t::cal_elemNum(vector<Bond_t>& _bonds){
	int rBondNum = 0;
	int cBondNum = 0;
	for(int b = 0; b < _bonds.size(); b++)
		if(_bonds[b].type == BD_ROW)
			rBondNum++;
		else if(_bonds[b].type == BD_COL)
			cBondNum++;
	Qnum_t qnum(0, 0);
	int dim;
	map<Qnum_t,int> row_QnumMdim;
	vector<int> row_offs(rBondNum, 0);
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

	map<Qnum_t,int> col_QnumMdim;
	vector<int> col_offs(cBondNum, 0);
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
	map<Qnum_t,int>::iterator it;
	map<Qnum_t,int>::iterator it2;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		_elemNum += it->second * it2->second;
	}
	return _elemNum;
}

Network_t::Network_t(): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
}

Network_t::Network_t(vector<SyTensor_t*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	for(int i = 0; i < tens.size(); i++)
		add(*tens[i]);
}

Network_t::Network_t(const string& fname): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	fromfile(fname);
	int Tnum = label_arr.size() - 1;
	order.assign(Tnum, 0);
	leafs.assign(Tnum, NULL);
	for(int i = 0; i < order.size(); i++)
		order[i] = i;
}

Network_t::Network_t(const string& fname, vector<SyTensor_t*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	fromfile(fname);
	assert((label_arr.size() - 1) == tens.size());
	for(int i = 0; i < tens.size(); i++){
		if(tens[i]->name.length() > 0){
			cout << tens[i]->name << ", "<< names[i]<<endl;
			assert(tens[i]->name == names[i]);
		}
		else
			tens[i]->setName(names[i]);
		tens[i]->addLabel(label_arr[i]);
		Node_t* ndp = new Node_t(*tens[i]);
		order.push_back(leafs.size());
		leafs.push_back(ndp);
	}
}

void Network_t::fromfile(const string& fname){
	string str;
	ifstream infile;
	infile.open (fname.c_str());
	int lnum = 0;
	int MAXLINES = 1000;
	int pos;
	int endpos;
	string tar("1234567890-");
	while(lnum < MAXLINES){
		getline(infile, str); // Saves the line in STRING.
		if(infile.eof())
			break;
		pos = str.find(":");
		assert(pos != string::npos);
		names.push_back(str.substr(0, pos));
		vector<int> labels;
		while((pos = str.find_first_of(tar, pos + 1)) != string::npos){
			endpos = str.find_first_not_of(tar, pos + 1);
			string label;
			if(endpos == string::npos)
				label = str.substr(pos);
			else
				label = str.substr(pos, endpos - pos);
			char* pEnd;
			labels.push_back(strtol(label.c_str(), &pEnd, 10));
			pos = endpos;
			if(pos == string::npos)
				break;
		}
		label_arr.push_back(labels);
		lnum ++;
	}
	assert(lnum < MAXLINES);	
	int numT = names.size() - 1;
	assert(names[numT] == "TOUT");
	infile.close();
}


Node_t* Network_t::add(SyTensor_t& SyTp){
	assert(label_arr.size() == 0);
	Node_t* ndp = new Node_t(SyTp);
	order.push_back(leafs.size());
	leafs.push_back(ndp);
	return ndp;
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
	load = false;
}

void Network_t::construct(){
	for(int i = 0; i < order.size(); i++)
		matching(leafs[order[i]], root);
	load = true;
}

//void Network_t::optimize(int num){
//	assert(false);//not a ready function
//}

SyTensor_t Network_t::launch(const string& _name){
	if(!load)
		construct();
	SyTensor_t SyT = merge(root);
	if(label_arr.size() > 0)
		SyT.reshape(label_arr[label_arr.size() - 1], 0);
	SyT.setName(_name);
	return SyT;
		
}

SyTensor_t Network_t::launch(int* outLabels, int Rnum, const string& _name){
	if(!load)
		construct();
	SyTensor_t SyT = merge(root);
	SyT.reshape(outLabels, Rnum);
	SyT.setName(_name);
	return SyT;
}

SyTensor_t Network_t::launch(vector<int>& outLabels, int Rnum, const string& _name){
	if(!load)
		construct();
	SyTensor_t SyT = merge(root);
	SyT.reshape(outLabels, Rnum);
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
	else
		if(nd->right->T == NULL){
			SyTensor_t rhtT = merge(nd->right);
			return *(nd->left->T) * rhtT;
		}
		else{
			return *(nd->left->T) * *(nd->right->T);
		}
}

Network_t::~Network_t(){
	if(load)
		destruct();
	for(int i = 0; i < leafs.size(); i++)
		delete leafs[i];
}

void Network_t::preprint(ostream& os, Node_t* nd, int layer){
	if(nd == NULL)
		return;
	for(int i = 0; i < layer; i++)
		os<<"|   ";
	if(nd->T)
		os<<nd->T->name << "(" << nd->elemNum << "): ";
	else
		os<<"*("<<nd->elemNum<<"): ";
	for(int i = 0; i < nd->labels.size(); i++)
		os<< nd->labels[i] << ", ";
	os<<endl;
	preprint(os, nd->left, layer+1);
	preprint(os, nd->right, layer+1);
}

ostream& operator<< (ostream& os, Network_t& net){
	net.preprint(os, net.root, 0);
	return os;
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
