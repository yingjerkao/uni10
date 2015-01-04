/****************************************************************************
*  @file Network.h
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
*  @brief Header file for Newtork class
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef NETWORK_H
#define NETWORK_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
#include <map>
#include <stdexcept>
#include <sstream>
//Bond property
#include <uni10/data-structure/uni10_struct.h>
#include <uni10/data-structure/Bond.h>
namespace uni10 {
class Network {
public:
    Network(const std::string& fname, const std::vector<UniTensor*>& tens);
    Network(const std::string& fname);
    ~Network();
    void putTensor(size_t idx, const UniTensor& UniT, bool force=true);
    void putTensor(size_t idx, const UniTensor* UniT, bool force=true); //if force is true, force replace without change the all network
    void putTensor(const std::string& name, const UniTensor& UniT, bool force=true);
    void putTensor(const std::string& name, const UniTensor* UniT, bool force=true);    //if force is true, force replace without change the all network
    void putTensorT(const std::string& nameT, const UniTensor& UniT, bool force=true);
    void putTensorT(const std::string& nameT, const UniTensor* UniT, bool force=true);
    UniTensor launch(const std::string& name="");
    std::string profile(bool print=true);
    friend std::ostream& operator<< (std::ostream& os, Network& net);
private:
    void preprint(std::ostream& os, Node* nd, int layer)const;  //pre-order print
    std::vector<std::string> names;
    std::map<std::string, size_t> name2pos;
    std::vector< std::vector<int> > label_arr;
    std::vector< int > Rnums;
    std::vector<Node*> leafs;
    std::vector<UniTensor*> tensors;
    std::vector< std::vector<_Swap> > swaps_arr;
    std::vector<bool> swapflags;
    std::vector<int> conOrder;  //contraction order;
    std::vector<int> order; //add order
    std::vector<int> brakets;   //add order
    Node* root;
    bool load;  //whether or not the network is ready for contraction, construct=> load=true, destruct=>load=false
    int times;  //construction times
    int tot_elem;   //total memory ussage
    int max_elem;   //maximum
    void construct();
    void destruct();
    void matching(Node* sbj, Node* tar);
    void branch(Node* sbj, Node* tar);
    UniTensor merge(Node* nd);
    void clean(Node* nd);
    void fromfile(const std::string& fname);
    void findConOrd(Node* nd);
    void addSwap();
    int rollcall();
    size_t sum_of_memory_usage();
    size_t max_tensor_elemNum();
    size_t memory_requirement();
    void _max_tensor_elemNum(Node* nd, size_t& max_num, Node& max_nd) const;
    size_t _sum_of_tensor_elem(Node* nd) const;
    size_t _elem_usage(Node* nd, size_t& usage, size_t& max_usage)const;
};
};  /* namespace uni10 */
#endif /* NETWORK_H */
