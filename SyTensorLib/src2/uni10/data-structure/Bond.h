#ifndef BOND_H
#define BOND_H

#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
#include <map>
#include <uni10/datatype.hpp>

namespace uni10{
enum bondType{
	BD_IN = 1,
	BD_OUT = -1
};
class UniTensor;
class Bond {
	public:
		Bond(bondType _type, int dim);
		Bond(bondType, const std::vector<Qnum>& qnums);
		Bond(const Bond& _b);
		void assign(bondType, int dim);
		void assign(bondType, const std::vector<Qnum>& qnums);
		bondType type()const;
		int dim()const;
		friend class UniTensor;
		friend class Node;
		friend std::ostream& operator<< (std::ostream& os, const Bond& b);
		friend bool operator== (const Bond& b1, const Bond& b2);
		void change(bondType tp);
		Bond& combine(const Bond bd);
		friend Bond combine(bondType tp, const std::vector<Bond>& bds);
		friend Bond combine(const std::vector<Bond>& bds);
		std::map<Qnum, int> degeneracy()const;
		std::vector<Qnum> Qlist()const;
		~Bond();
	private:
		void setting(const std::vector<Qnum>& qnums);
		bondType m_type;
		int m_dim;
		std::vector<Qnum>Qnums;	//Quantum numbers
		std::vector<int>Qdegs;	//Degeneracy in each quantum sector
		std::vector<int>offsets;	
};
Bond combine(bondType tp, const std::vector<Bond>& bds);
Bond combine(const std::vector<Bond>& bds);
};
#endif /* BOND_H */
