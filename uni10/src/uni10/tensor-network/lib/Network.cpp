/****************************************************************************
*  @file Network.cpp
*  @license
*    Universal Tensor Network Library
*    Copyright (c) 2013-2014
*    National Taiwan University
*    National Tsing-Hua University

*
*    This file is part of Uni10, the Universal Tensor Network Library.
*
*    Uni10 is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Lesser General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Uni10 is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public License
*    along with Uni10.  If not, see <http://www.gnu.org/licenses/>.
*  @endlicense
*  @brief Implementation file for Node and Network classes 
*  @author Yun-Da Hsieh, Ying-Jer Kao
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#include <algorithm>
#include <uni10/tools/uni10_tools.h>
#include <uni10/tensor-network/UniTensor.h>
#include <uni10/tensor-network/Network.h>


namespace uni10{
Node::Node(): T(NULL), elemNum(0), parent(NULL), left(NULL), right(NULL), point(0){
}

Node::Node(UniTensor* Tp): T(Tp), elemNum(Tp->m_elemNum), labels(Tp->labels), bonds(Tp->bonds), name(Tp->name), parent(NULL), left(NULL), right(NULL), point(0){
  if(!(Tp->status & Tp->HAVEBOND)){
    std::ostringstream err;
    err<<"Cannot create node of a network from tensor without bond.";
    throw std::runtime_error(exception_msg(err.str()));
  }
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
	size_t dim;
	std::map<Qnum,size_t> row_QnumMdim;
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

	std::map<Qnum,size_t> col_QnumMdim;
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
	size_t _elemNum = 0;
	std::map<Qnum,size_t>::iterator it;
	std::map<Qnum,size_t>::iterator it2;
	for ( it2 = col_QnumMdim.begin() ; it2 != col_QnumMdim.end(); it2++ ){
		it = row_QnumMdim.find(it2->first);
		_elemNum += it->second * it2->second;
	}
	return _elemNum;
}


Network::Network(const std::string& fname): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
  try{
    fromfile(fname);
    int Tnum = label_arr.size() - 1;
    swapflags.assign(Tnum, false);
    std::vector<_Swap> swaps;
    swaps_arr.assign(Tnum, swaps);
    leafs.assign(Tnum, (Node*)NULL);
    tensors.assign(Tnum, (UniTensor*)NULL);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Network::Network(std::string&):");
  }
}

Network::Network(const std::string& fname, const std::vector<UniTensor*>& tens): root(NULL), load(false), times(0), tot_elem(0), max_elem(0){
  try{
    fromfile(fname);
    if(!((label_arr.size() - 1) == tens.size())){
      std::ostringstream err;
      err<<"The size of the input vector of tensors does not match for the number of tensors in the network file '"<<fname<<"'.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    int Tnum = tens.size();
    swapflags.assign(Tnum, false);
    std::vector<_Swap> swaps;
    swaps_arr.assign(Tnum, swaps);
    leafs.assign(Tnum, (Node*)NULL);
    tensors.assign(Tnum, (UniTensor*)NULL);
    for(int i = 0; i < Tnum; i++){
      if(!(tens[i]->RBondNum == Rnums[i])){
        std::ostringstream err;
        err<<"The number of in-coming bonds does not match with the tensor '"<<names[i]<<"' specified in network file '"<<fname<<"'.";
        throw std::runtime_error(exception_msg(err.str()));
      }
      UniTensor* ten = new UniTensor(*(tens[i]));
      ten->setName(names[i]);
      ten->setLabel(label_arr[i]);
      tensors[i] = ten;
      Node* ndp = new Node(ten);
      leafs[i] = ndp;
    }
    construct();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In constructor Network::Network(std::string&, std::vector<UniTensor*>&):");
  }
}

void Network::fromfile(const std::string& fname){//names, name2pos, label_arr, Rnums, order, brakets
	std::string str;
	std::ifstream infile;
	infile.open (fname.c_str());
	if(!(infile.is_open())){
      std::ostringstream err;
      err<<"Error in opening file '" << fname <<"'.";
      throw std::runtime_error(exception_msg(err.str()));
  }
	int lnum = 0;
	int MAXLINES = 1000;
	int pos = 0;
	int endpos = 0;
	std::string tar("1234567890-");
	std::vector<std::string> ord;
	while(lnum < MAXLINES){
		getline(infile, str); // Saves the line in STRING.
    endpos = 0;
		if(infile.eof())
			break;
		pos = str.find(":");
		if(pos == std::string::npos)
			break;
		std::string name = str.substr(0, pos);
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
        if(!(sum == 0)){
          std::ostringstream err;
          err<<"Error in the network file '"<<fname<<"'. There are imbalance brackets when specifying the contraction order.";
          throw std::runtime_error(exception_msg(err.str()));
        }
			}
			break;
		}
		name2pos[name] = names.size();
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
	std::vector<bool> found(numT, false);
  if(!(names[numT] == "TOUT")){
    std::ostringstream err;
    err<<"Error in the network file '"<<fname<<"'. Missing TOUT tensor. One must specify the bond labels for the resulting tensor by giving TOUT tag.\n  Hint: If the resulting tensor is a scalar(0-bond tensor), leave the TOUT empty like 'TOUT: '";
    throw std::runtime_error(exception_msg(err.str()));
  }
	if(!(names.size() > 2)){
    std::ostringstream err;
    err<<"Error in the network file '"<<fname<<"'. There must be at least two tensors in a tensor network.";
    throw std::runtime_error(exception_msg(err.str()));
  }
	order.assign(numT, 0);
	if(ord.size() > 0){
    if(!(ord.size() == numT)){
      std::ostringstream err;
      err<<"Error in the network file '"<<fname<<"'. Some tensors are missing in the contraction order.";
      throw std::runtime_error(exception_msg(err.str()));
    }
		std::map<std::string, size_t>::iterator it;
		for(int i = 0; i < numT; i++){
			it = name2pos.find(ord[i]);
      if(!(it != name2pos.end())){
        std::ostringstream err;
        err<<"Error in the network file '"<<fname<<"'. '"<<ord[i]<<"' in the contraction order is not in the list of tensors above.";
        throw std::runtime_error(exception_msg(err.str()));
      }
			order[i] = it->second;
			if(!(found[order[i]] == false)){
        std::ostringstream err;
        err<<"Error in the network file '"<<fname<<"'. '"<<ord[i]<<"' appears more than once in the contraction order.";
        throw std::runtime_error(exception_msg(err.str()));
      }
			found[order[i]] = true;
		}
	}
	else{
		for(int i = 0; i < numT; i++)
			order[i] = i;
	}
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
        if(leafs[order[cnt]] == NULL){
          std::ostringstream err;
          err<<"(((Tensor '"<<names[order[cnt]]<<"' has not yet been given.\n  Hint: Use addTensor() to add a tensor to a network.\n";
          throw std::runtime_error(exception_msg(err.str()));
        }
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
			if(leafs[order[i]] == NULL){
          std::ostringstream err;
          err<<"Tensor '"<<names[order[i]]<<"' has not yet been given.\n  Hint: Use putTensor() to add a tensor to a network.\n";
          throw std::runtime_error(exception_msg(err.str()));
      }
			matching(leafs[order[i]], root);
		}
	}
  int Tnum = label_arr.size() - 1;
  if(root->labels.size() == label_arr[Tnum].size()){
    for(int l = 0; l < root->labels.size(); l++){
      bool found = false;
      for(int t = 0; t < label_arr[Tnum].size(); t++)
        if(root->labels[l] == label_arr[Tnum][t]){
          found = true;
          break;
        }
      if(!found){
        std::ostringstream err;
        err<<"Error when constructing the network. The labels of the resulting tensor, ( ";
        for(int i = 0; i < root->labels.size(); i++)
          err<<root->labels[i]<<" ";
        err<<"), do not match with the labels of 'TOUT' in the network file";
        throw std::runtime_error(exception_msg(err.str()));
      }
    }

  }
  else{
    std::ostringstream err;
    err<<"Error when constructing the network. The bond number of the resulting tensor is different from the bond number of 'TOUT'";
    throw std::runtime_error(exception_msg(err.str()));
  }
	addSwap();
	load = true;
}

void Network::putTensor(size_t idx, const UniTensor* UniT, bool force){
  try{
    if(!(idx < (label_arr.size()-1))){
      std::ostringstream err;
      err<<"Index exceeds the number of the tensors in the list of network file.";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if((!force) && load){
      destruct();
    }
    if(!(UniT->RBondNum == Rnums[idx])){
      std::ostringstream err;
      err<<"The number of in-coming bonds does not match with the tensor '"<<names[idx]<<"' specified in network file";
      throw std::runtime_error(exception_msg(err.str()));
    }
    if(leafs[idx] != NULL){
      *(tensors[idx]) = *UniT;
      tensors[idx]->setLabel(label_arr[idx]);
      tensors[idx]->setName(names[idx]);
      swapflags[idx] = false;
    }
    else{
      UniTensor* ten = new UniTensor(*UniT);
      ten->setName(names[idx]);
      ten->setLabel(label_arr[idx]);
      tensors[idx] = ten;
      Node* ndp = new Node(ten);
      leafs[idx] = ndp;
    }
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensor(size_t, uni10::UniTensor*, bool=true):");
  }
}

void Network::putTensor(size_t idx, const UniTensor& UniT, bool force){
  try{
    putTensor(idx, &UniT, force);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensor(size_t, uni10::UniTensor&, bool=true):");
  }
}

void Network::putTensor(const std::string& name, const UniTensor* UniT, bool force){
  try{
    std::map<std::string, size_t>::const_iterator it = name2pos.find(name);
    if(!(it != name2pos.end())){
      std::ostringstream err;
      err<<"There is no tensor named '"<<name<<"' in the network file";
      throw std::runtime_error(exception_msg(err.str()));
    }
    putTensor(it->second, UniT, force);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensor(std::string&, uni10::UniTensor*, bool=true):");
  }
}

void Network::putTensor(const std::string& name, const UniTensor& UniT, bool force){
  try{
    putTensor(name, &UniT, force);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensor(std::string&, uni10::UniTensor&, bool=true):");
  }
}

void Network::putTensorT(const std::string& nameT, const UniTensor* UniT, bool force){
  try{
    std::map<std::string, size_t>::const_iterator itT = name2pos.find(nameT);
    if(!(itT != name2pos.end())){
      std::ostringstream err;
      err<<"There is no tensor named '"<<nameT<<"' in the network file";
      throw std::runtime_error(exception_msg(err.str()));
    }
    UniTensor transT = *UniT;
    transT.transpose();
    putTensor(itT->second, &transT, force);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensorT(std::string&, uni10::UniTensor*, bool=true):");
  }
}

void Network::putTensorT(const std::string& nameT, const UniTensor& UniT, bool force){
  try{
    putTensorT(nameT, &UniT, force);
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::putTensorT(std::string&, uni10::UniTensor&, bool=true):");
  }
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
      if(!(tar->left != NULL && tar->right != NULL)){
	std::ostringstream err;
	err<<"Fatal error(code = N1). Please contact the developer of the uni10 library.";
	throw std::runtime_error(exception_msg(err.str()));
      }
      float tar_p = tar->point;
      float lft_p = 0, rht_p = 0;
      if((lft_p = sbj->metric(tar->left)) > tar_p || (rht_p = sbj->metric(tar->right)) > tar_p){	//go deeper comparison to the children
	if(lft_p > rht_p)
	  matching(sbj, tar->left);
	else
	  matching(sbj, tar->right);
      }
      else	//contract
	branch(sbj, tar);
    }
    else	//contract
      branch(sbj, tar);
  }
  else{	//contract
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
  try{
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
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::launch(std::string&):");
    return UniTensor();
  }
}

UniTensor Network::merge(Node* nd){
  if(nd->left->T == NULL){
    UniTensor lftT = merge(nd->left);
    if(nd->right->T == NULL){
      UniTensor rhtT = merge(nd->right);
      return contract(lftT, rhtT, true);
    }
    else{
      return contract(lftT, *(nd->right->T), true);
    }
  }
  else{
    if(nd->right->T == NULL){
      UniTensor rhtT = merge(nd->right);
      return contract(*(nd->left->T), rhtT, true);
    }
    else{
      return contract(*(nd->left->T), *(nd->right->T), true);
    }
  }
}

Network::~Network(){
  try{
    if(load)
      destruct();
    for(int i = 0; i < leafs.size(); i++)
      delete leafs[i];
    for(int i = 0; i < tensors.size(); i++)
      delete tensors[i];
  }
  catch(const std::exception& e){
    propogate_exception(e, "In destructor Network::~Network():");
  }
}

int Network::rollcall(){
  if(!load){
    for(int i = 0; i < leafs.size(); i++)
      if(leafs[i] == NULL){
        return i;
      }
    construct();
  }
  return -1;
}

size_t Network::sum_of_memory_usage(){
  if(rollcall() >= 0)
    return 0;
  return _sum_of_tensor_elem(root) * sizeof(Real);
}

size_t Network::_sum_of_tensor_elem(Node* nd) const{
  if(nd == NULL)
    return 0;
  return nd->elemNum + _sum_of_tensor_elem(nd->left) + _sum_of_tensor_elem(nd->right);
}

size_t Network::memory_requirement(){
  if(rollcall() >= 0)
    return 0;
  size_t usage = 0;
  for(int i = 0; i < leafs.size(); i++)
    usage += leafs[i]->elemNum;
  usage *= 2;
  size_t max_usage = 0;
  _elem_usage(root, usage, max_usage);
  return max_usage * sizeof(Real);
}

size_t Network::_elem_usage(Node* nd, size_t& usage, size_t& max_usage)const{
  if(nd == NULL)
    return 0;
  size_t child_usage = _elem_usage(nd->left, usage, max_usage) + _elem_usage(nd->right, usage, max_usage);
  usage += nd->elemNum;
  max_usage = usage > max_usage ? usage : max_usage;
  usage -= child_usage;
  return nd->elemNum;
}

size_t Network::max_tensor_elemNum(){
  if(rollcall() >= 0)
    return 0;
  size_t max_num = 0;
  Node max_nd;
  _max_tensor_elemNum(root, max_num, max_nd);
  return max_num;
}

void Network::_max_tensor_elemNum(Node* nd, size_t& max_num, Node& max_nd) const{
  if(nd == NULL)
    return;
  _max_tensor_elemNum(nd->left, max_num, max_nd);
  _max_tensor_elemNum(nd->right, max_num, max_nd);
  if(nd->elemNum > max_num){
    max_num = nd->elemNum;
    max_nd = *nd;
  }
}

std::string Network::profile(bool print){
  try{
    std::ostringstream os;
    int miss;
    if((miss = rollcall()) >= 0){
      os<<"\nTensor '"<<names[miss]<<"' has not yet been given!\n\n";
      if(print){
        std::cout<<os.str();
        return "";
      }
      return os.str();
    }
    os<<"\n===== Network profile =====\n";
    os<<"Memory Requirement: "<<memory_requirement()<<std::endl;
    //os<<"Sum of memory usage: "<<sum_of_memory_usage()<<std::endl;
    size_t max_num = 0;
    Node max_nd;
    _max_tensor_elemNum(root, max_num, max_nd);
    os<<"Maximun tensor: \n";
    os<<"  elemNum: "<<max_num<<"\n  "<<max_nd.labels.size()<<" bonds and labels: ";
    for(int i = 0; i < max_nd.labels.size(); i++)
      os<< max_nd.labels[i] << ", ";
    os<<std::endl;
    os<<"===========================\n\n";
    if(print){
      std::cout<<os.str();
      return "";
    }
    return os.str();
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function Network::profile():");
    return "";
  }
}

void Network::preprint(std::ostream& os, Node* nd, int layer)const{
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
  try{
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
    if(net.rollcall() < 0)
      net.preprint(os, net.root, 0);
    else
      os<<"\nSome tensors have not yet been given!\n\n";
  }
  catch(const std::exception& e){
    propogate_exception(e, "In function operator<<(std::ostream&, uni10::Network&):");
  }
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
    if(!found){
      std::ostringstream err;
      err<<"Fatal error(code = N2). Please contact the developer of the uni10 library.";
      throw std::runtime_error(exception_msg(err.str()));
    }
	}
	findConOrd(nd->left);
	findConOrd(nd->right);
}

void Network::addSwap(){
	int Tnum = leafs.size();
	findConOrd(root);
	if(!(Tnum == conOrder.size())){
      std::ostringstream err;
      err<<"Fatal error(code = N3). Please contact the developer of the uni10 library.";
      throw std::runtime_error(exception_msg(err.str()));
  }
	//int tenOrder[conOrder.size()];
  std::vector<int> tenOrder = conOrder;
	//memcpy(tenOrder, &(conOrder[0]), Tnum * sizeof(int));
	std::vector<_Swap> tenSwaps = recSwap(tenOrder);
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
