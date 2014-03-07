#ifndef UNI10_STRUCT_H
#define UNI10_STRUCT_H
#include <string>
#include <stdint.h>
namespace uni10{
typedef struct{
	int b1; 
	int b2; 
}_Swap;
class UniTensor;
class Bond;
class Node{
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
		UniTensor* T;	//if T != NULL, it is leaf node
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
