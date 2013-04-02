#include "Block.h"
ostream& operator<< (ostream& os, const Block_t& b){
	os << ": " << b.Rnum << " x " << b.Cnum << " = " << b.Rnum * b.Cnum;
	return os;
}
int Block_t::row(){
	return Rnum;
}
int Block_t::col(){
	return Cnum;
}
bool operator== (const Block_t& b1, const Block_t& b2){
	return (b1.qnum == b2.qnum) && (b1.Rnum == b2.Rnum) && (b1.Cnum == b2.Cnum);
}

