#ifndef UNI10_STRUCT_H
#define UNI10_STRUCT_H
typedef struct{
	int b1; 
	int b2; 
}_Swap;
class SyTensor_t;
class Bond_t;
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

#endif /* UNI10_STRUCT_H */
