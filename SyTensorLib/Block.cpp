#include "Block.h"
ostream& operator<< (ostream& os, const Block_t& b){
	os << "--- " << b.qnum<< ": " << b.Rnum << " x " << b.Cnum << " = " << b.Rnum * b.Cnum << " ---\n\n";
	for(int r = 0; r < b.Rnum; r++){
		for(int c = 0; c < b.Cnum; c++)
			cout<< setw(7) << setprecision(3) << b.elem[r * b.Cnum + c];
		os << "\n\n";
	}
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

