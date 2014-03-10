#include <algorithm>
#include <uni10/tensor-network/Network.h>
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/tools/uni10_tools.h>
namespace uni10{
Node::Node(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL), point(0){
}

Node::Node(UniTensor* Tp): T(Tp), elemNum(Tp->m_elemNum), labels(Tp->labels), bonds(Tp->bonds), name(Tp->name), parent(NULL), left(NULL), right(NULL), point(0){	
	assert(Tp->status & Tp->HAVEBOND);
	//assert(Tp->status & Tp->HAVELABEL);
}

Node::Node(const Node& nd): T(nd.T), elemNum(nd.elemNum), labels(nd.labels), bonds(nd.bonds), parent(nd.parent), left(nd.left), right(nd.right), point(nd.point){	
}

Node::Node(std::vector<Bond>& _bonds, std::vector<int>& _labels): T(NULL), labels(_labels), bonds(_bonds), parent(NULL), left(NULL), right(NULL), point(0){	
	elemNum = cal_elemNum(bonds);
}

Node::~Node(){
}

void Node::delink(){
	parent = NULL;
	left = NULL;
	right = NULL;
	point = 0;
}

Node Node::contract(Node* nd){
	int AbondNum = bonds.size();
	int BbondNum = nd->bonds.size();
	std::vector<Bond> cBonds;
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
		cBonds[a].change(BD_IN);
	for(int a = 0; a < cBondNum; a++)
		cBonds[rBondNum + a].change(BD_OUT);

	//Node par(cBonds, newLabelC);
	return Node(cBonds, newLabelC);
}

float Node::metric(Node* nd){	//Bigger is better
	int AbondNum = bonds.size();
	int BbondNum = nd->bonds.size();
	std::vector<Bond> cBonds;
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
		cBonds[a].change(BD_IN);
	for(int a = 0; a < cBondNum; a++)
		cBonds[rBondNum + a].change(BD_OUT);
	int64_t newElemNum = cal_elemNum(cBonds);
	return float(elemNum + nd->elemNum) / newElemNum;
}

int64_t Node::cal_elemNum(std::vector<Bond>& _bonds){
	int rBondNum = 0;
	int cBondNum = 0;
	for(int b = 0; b < _bonds.size(); b++)
		if(_bonds[b].type() == BD_IN)
			rBondNum++;
		else if(_bonds[b].type() == BD_OUT)
			cBondNum++;
	Qnum qnum(0, PRT_EVEN);
	int dim;
	std::map<Qnum,int> row_QnumMdim;
	std::vector<int> row_offs(rBondNum, 0);
	if(rBondNum){
		while(1){
			qnum.assign();
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
		qnum.assign();
		row_QnumMdim[qnum] = 1;
	}

	std::map<Qnum,int> col_QnumMdim;
	std::vector<int> col_offs(cBondNum, 0);
	if(cBondNum){
		while(1){
			qnum.assign();
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
		qnum.assign();
		if(row_QnumMdim.find(qnum) != row_QnumMdim.end())
			col_QnumMdim[qnum] = 1;
	}
	int64_t _elemNum = 0;
	std::map<Qnum,int>::iterator it;
	std::map<Qnum,int>::iterator it2;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		_elemNum += it->second * it2->second;
	}
	return _elemNum;
}

/*
Network::Network(): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
}
*/
/*
Network::Network(std::vector<UniTensor*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	for(int i = 0; i < tens.size(); i++)
		add(tens[i]);
}
*/

Network::Network(const std::string& fname): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
	fromfile(fname);
	int Tnum = label_arr.size() - 1;
	swapflags.assign(Tnum, false);
	std::vector<_Swap> swaps;
	swaps_arr.assign(Tnum, swaps);
	leafs.assign(Tnum, NULL);
	tensors.assign(Tnum, NULL);
}

Network::Network(const std::string& fname, const std::vector<UniTensor*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
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
		UniTensor* ten = new UniTensor(*(tens[i]));
		ten->setName(names[i]);
		ten->addLabel(label_arr[i]);
		tensors[i] = ten;
		Node* ndp = new Node(ten);
		leafs[i] = ndp;
	}
	construct();
}

void Network::fromfile(const std::string& fname){//names, label_arr, Rnums, order, brakets
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
		//boost::algorithm::trim(name);
		trim(name);
		if(name == "ORDER"){
			std::string bra("(");
			if(str.find(bra, pos+1) == std::string::npos){
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
			}
			else{
				std::string del(" ,;*()");
				endpos = str.find_first_of(del, pos + 1);
				while(((pos = str.find_first_not_of(del, pos + 1)) != std::string::npos)){
					std::string bras = str.substr(endpos, pos - endpos);
					int minus = std::count(bras.begin(), bras.end(), ')');
					int plus = std::count(bras.begin(), bras.end(), '(');
					if(minus)
						brakets.push_back(-minus);
					if(plus)
						brakets.push_back(plus);
					endpos = str.find_first_of(del, pos + 1);
					if(endpos == std::string::npos){
						ord.push_back(str.substr(pos));
						brakets.push_back(0);
					}
					else{
						ord.push_back(str.substr(pos, endpos - pos));
						brakets.push_back(0);
					}
					pos = endpos;
					if(pos == std::string::npos)
						break;
				}
				if(endpos != std::string::npos){
					std::string bras = str.substr(endpos);
					int minus = std::count(bras.begin(), bras.end(), ')');
					if(minus)
						brakets.push_back(-minus);
				}
				int sum = 0;
				for(int i = 0; i < brakets.size(); i++){
					sum += brakets[i];
				}
				assert(sum == 0);
			}
			break;
		}
		names.push_back(name);
		std::vector<int> labels;
		int Rnum = -1;
		int cnt = 0;
		int tmp;
		while((pos = str.find_first_of(tar, pos + 1)) != std::string::npos){
			if(Rnum == -1){
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
			if(Rnum == -1)
				cnt++;
			if(pos == std::string::npos)
				break;
		}
		label_arr.push_back(labels);
		if(Rnum == -1)
			Rnum = labels.size();
		Rnums.push_back(Rnum);
		lnum ++;
	}
	int numT = names.size() - 1;
	assert(names[numT] == "TOUT");
	assert(names.size() > 0);
	order.assign(numT, 0);
	if(ord.size() > 0){
		assert(ord.size() == numT);
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
	/*
	std::cout<<"Brakets: ";
	for(int i = 0; i < brakets.size(); i++)
		std::cout<<brakets[i]<<", ";
	std::cout<<"\norder: ";
	for(int i = 0; i < numT; i++)
		std::cout<<order[i]<<", ";
	std::cout<<std::endl;
	*/
	infile.close();
}

void Network::construct(){
	if(brakets.size()){
		std::vector<Node*> stack(leafs.size(), NULL);
		int cursor = 0;
		int cnt = 0;
		for(int i = 0; i < brakets.size(); i++){
			if(brakets[i] < 0){
				for(int c = 0; c < -brakets[i]; c++){
					Node* par = new Node(stack[cursor - 2]->contract(stack[cursor - 1]));
					par->left = stack[cursor - 2];
					par->right = stack[cursor - 1];
					par->left->parent = par;
					par->right->parent = par;
					stack[cursor - 2] = par;
					cursor--;
					if(cursor < 2)	//prevent breakdown because of extra braket
						break;
				}
			}
			else if(brakets[i] == 0){
				stack[cursor] = leafs[order[cnt]];
				cnt++;
				cursor++;
			}
		}
		while(cursor > 1){//for imcomplete brakets
			Node* par = new Node(stack[cursor - 2]->contract(stack[cursor - 1]));
			par->left = stack[cursor - 2];
			par->right = stack[cursor - 1];
			par->left->parent = par;
			par->right->parent = par;
			stack[cursor - 2] = par;
			cursor--;
		}
		root = stack[0];
	}
	else{
		for(int i = 0; i < order.size(); i++){
			assert(leafs[order[i]] != NULL);
			matching(leafs[order[i]], root);
		}
	}
	addSwap();
	load = true;
}

Node* Network::putTensor(int idx, const UniTensor* UniT, bool force){
	assert(label_arr.size() > 0 && idx >= 0 && idx < (label_arr.size()-1));
	if((!force) && load)
		destruct();
	if(UniT->name.length() > 0)
		assert(UniT->name == names[idx]);
	assert(UniT->RBondNum == Rnums[idx]);

	if(leafs[idx] != NULL){
		assert(tensors[idx]->similar(*UniT));
		*(tensors[idx]) = *UniT;
		tensors[idx]->addLabel(label_arr[idx]);
		tensors[idx]->setName(names[idx]);
		swapflags[idx] = false;
	}
	else{
		UniTensor* ten = new UniTensor(*UniT);
		ten->setName(names[idx]);
		ten->addLabel(label_arr[idx]);
		tensors[idx] = ten;
		Node* ndp = new Node(ten);
		leafs[idx] = ndp;
	}
	return leafs[idx];
}

void Network::branch(Node* sbj, Node* tar){
	Node* par = new Node(tar->contract(sbj));
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

void Network::matching(Node* sbj, Node* tar){
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

void Network::clean(Node* nd){
	if(nd->T != NULL)	//leaf
		return;
	clean(nd->left);
	clean(nd->right);
	delete nd;
}

void Network::destruct(){
	clean(root);
	root = NULL;
	for(int i = 0; i < leafs.size(); i++)
		leafs[i]->delink();
	conOrder.clear();
	for(int t = 0; t < tensors.size(); t++){
		if(Qnum::isFermionic() && swapflags[t]){
			tensors[t]->addGate(swaps_arr[t]);
			swapflags[t] = false;
		}
		swaps_arr[t].clear();
	}
	load = false;
}


UniTensor Network::launch(const std::string& _name){
	if(!load)
		construct();
	for(int t = 0; t < tensors.size(); t++)
		if(Qnum::isFermionic() && !swapflags[t]){
			tensors[t]->addGate(swaps_arr[t]);
			swapflags[t] = true;
		}
	UniTensor UniT = merge(root);
	int idx = label_arr.size() - 1;
	if(label_arr.size() > 0 && label_arr[idx].size() > 1)
		UniT.permute(label_arr[idx], Rnums[idx]);
	UniT.setName(_name);
	return UniT;
}

UniTensor Network::merge(Node* nd){
	if(nd->left->T == NULL){
		UniTensor lftT = merge(nd->left);
		if(nd->right->T == NULL){
			UniTensor rhtT = merge(nd->right);
			return lftT * rhtT;
		}
		else{
			return lftT * *(nd->right->T);
		}
	}
	else{
		if(nd->right->T == NULL){
			UniTensor rhtT = merge(nd->right);
			return *(nd->left->T) * rhtT;
		}
		else{
			return *(nd->left->T) * *(nd->right->T);
		}
	}
}

Network::~Network(){
	if(load)
		destruct();
	for(int i = 0; i < leafs.size(); i++)
		delete leafs[i];
	for(int i = 0; i < tensors.size(); i++)
		delete tensors[i];
}

void Network::preprint(std::ostream& os, Node* nd, int layer){
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

std::ostream& operator<< (std::ostream& os, Network& net){
	os<<std::endl;
	for(int i = 0; i < net.names.size(); i++){
		os<<net.names[i]<< ": ";
		if(net.Rnums[i])
			os<<"i[";
		for(int l = 0; l < net.Rnums[i]; l++){
			os<<net.label_arr[i][l];
			if(l < net.Rnums[i] - 1)
				os<<", ";
		}
		if(net.Rnums[i])
			os<<"] ";
		if(net.label_arr[i].size() - net.Rnums[i])
			os<<"o[";
		for(int l = net.Rnums[i]; l < net.label_arr[i].size(); l++){
			os<<net.label_arr[i][l];
			if(l < net.label_arr[i].size() - 1)
				os<<", ";
		}
		if(net.label_arr[i].size() - net.Rnums[i])
			os<<"]";
		os<<std::endl;
	}
	os<<std::endl;
	if(net.load)
		net.preprint(os, net.root, 0);

	return os;
}
std::ostream& operator<< (std::ostream& os, const Node& nd){
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

void Network::findConOrd(Node* nd){
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

void Network::addSwap(){
	int Tnum = leafs.size();
	findConOrd(root);
	assert(Tnum == conOrder.size());
	int tenOrder[conOrder.size()];
	memcpy(tenOrder, &(conOrder[0]), Tnum * sizeof(int));
	std::vector<_Swap> tenSwaps = recSwap(tenOrder, Tnum);
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
}; /* namespace uni10 */
