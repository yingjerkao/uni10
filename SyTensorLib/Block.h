#include <iostream>
#include <iomanip>
#include <assert.h>
using namespace std;
class SyTensor_t;
class Block_t {
	public:
		Block_t(): Rnum(0), Cnum(0), offset(0){
			//cout<<"Constructing Block...\n";
		}
		Block_t(int _Rnum, int _Cnum, int _offset): Rnum(_Rnum), Cnum(_Cnum), offset(_offset){
			//cout<<"Constructing Block...\n";
		}
		Block_t(const Block_t& _b): Rnum(_b.Rnum), Cnum(_b.Cnum), offset(_b.offset){
			//cout<<"Copying Block...\n";
		}
		~Block_t(){
			//cout<<"Destructing Block...\n";
		};
		friend class SyTensor_t;
		friend ostream& operator<< (ostream& os, const Block_t& b);
		friend ostream& operator<< (ostream& os, SyTensor_t& SyT);
	private:
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		int64_t offset;	//index of the first element of a block element in Tensor
};
