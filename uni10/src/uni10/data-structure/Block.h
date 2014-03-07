#ifndef BLOCK_H
#define BLOCK_H
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <stdint.h>
#include <uni10/datatype.hpp>

namespace uni10{
class UniTensor;
class Block{
	public:
		Block();
		Block(const Block& _b);
		~Block();
		friend class UniTensor;
		friend std::ostream& operator<< (std::ostream& os, const Block& b);
		friend std::ostream& operator<< (std::ostream& os, UniTensor& UniT);
		friend UniTensor operator* (UniTensor& Ta, UniTensor& Tb);
		friend bool operator== (const Block& b1, const Block& b2);
	private:
		Qnum qnum;
		double* elem;
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		int64_t offset;	//index of the first element of a block element in Tensor
};
};
#endif /* BLOCK_H */

