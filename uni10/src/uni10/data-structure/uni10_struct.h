/****************************************************************************
*  @file uni10_struct.h
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
*  @brief Generic header file for uni10 data strucutres
*  @author Yun-Da Hsieh
*  @date 2014-05-06
*  @since 0.1.0
*
*****************************************************************************/
#ifndef UNI10_STRUCT_H
#define UNI10_STRUCT_H
#include <cstdint>
#include <string>
#include <stdexcept>
#include <vector>
#include <uni10/data-structure/Block.h>
namespace uni10{
typedef struct{
	int b1;
	int b2;
}_Swap;
class UniTensor;
class Bond;
class Node {
public:
    Node();
    Node(UniTensor* Tp);
    Node(const Node& nd);
    Node(std::vector<Bond>& _bonds, std::vector<int>& _labels);
    ~Node();
    Node contract(Node* nd);
    float metric(Node* nd);
    friend std::ostream& operator<< (std::ostream& os, const Node& nd);
    friend class Network;
private:
    UniTensor* T;   //if T != NULL, it is leaf node
    std::vector<int> labels;
    std::vector<Bond> bonds;
    int64_t elemNum;
    std::string name;
    Node* parent;
    Node* left;
    Node* right;
    float point;
    int64_t cal_elemNum(std::vector<Bond>& _bonds);
    void delink();
};

}; /* namespace uni10 */

#endif /* UNI10_STRUCT_H */
