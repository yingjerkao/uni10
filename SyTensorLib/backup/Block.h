#include <iostream>
#include <iomanip>
#include <assert.h>
using namespace std;

class Block_t {
	public:
		Block_t(): Rnum(0), Cnum(0), offset(0){
			cout<<"Constructing Qnum...\n";
		}
		Block_t(int _Rnum, int _Cnum, int _offset): Rnum(_Rnum), Cnum(_Cnum), offset(_offset){
			cout<<"Constructing Qnum...\n";
		}
		~Block_t(){
			cout<<"Destructing Qnum...\n";
		};

		friend ostream& operator<< (ostream& os, const Block_t& b);
	private:
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		int64_t offset;	//index of the first element of a block element in Tensor
};

ostream& operator<< (ostream& os, const Block_t& b){
	os << "Block: " << b.Rnum << " x " << b.Cnum << " = " << b.Rnum * b.Cnum;
	return os;
}
