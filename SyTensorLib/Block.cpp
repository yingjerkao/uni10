#include "Block.h"
ostream& operator<< (ostream& os, const Block_t& b){
	os << ": " << b.Rnum << " x " << b.Cnum << " = " << b.Rnum * b.Cnum;
	return os;
}

