#ifndef BOND_H
#define BOND_H

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
#include <uni10/datatype.hpp>

//namespace uni10{
//namespace datatype{
enum bondType{
	BD_ROW = 1,
	BD_COL = -1
};
class SyTensor_t;
class Bond_t {
	public:
		Bond_t(bondType, std::vector<Qnum>& qnums);
		Bond_t(const Bond_t& _b):type(_b.type), dim(_b.dim), Qnums(_b.Qnums), Qdegs(_b.Qdegs), offsets(_b.offsets){
			//cout<<"Copying Bond "<< this <<" from " << &_b << endl;
		}
		void assign(bondType, std::vector<Qnum>& qnums);
		friend class SyTensor_t;
		friend class Node_t;
		friend std::ostream& operator<< (std::ostream& os, const Bond_t& b);
		friend std::ostream& operator<< (std::ostream& os, SyTensor_t& SyT);
		friend bool operator== (const Bond_t& b1, const Bond_t& b2);
		void change(bondType tp);
		void combine(Bond_t bd);
		friend Bond_t combine(bondType tp, const std::vector<Bond_t>& bds);
		friend Bond_t combine(const std::vector<Bond_t>& bds);
		~Bond_t();
	private:
		void setting(std::vector<Qnum>& qnums);
		bondType type;
		int dim;
		std::vector<Qnum>Qnums;	//Quantum numbers
		std::vector<int>Qdegs;	//Degeneracy in each quantum sector
		std::vector<int>offsets;	
};
//};
//};
#endif /* BOND_H */
