#include <iostream>
#include <iomanip>
#include <assert.h>
#include "Qnum.h"
using namespace std;
class SyTensor_t;
class Block_t {
	public:
		Block_t(): Rnum(0), Cnum(0), offset(0), elem(NULL){
			//cout<<"Constructing Block...\n";
		}
		Block_t(const Block_t& _b): qnum(_b.qnum), Rnum(_b.Rnum), Cnum(_b.Cnum), offset(_b.offset), elem(_b.elem){
			//cout<<"Copying Block...\n";
		}
		~Block_t(){
			//cout<<"Destructing Block...\n";
		};
		int row();
		int col();
		friend class SyTensor_t;
		friend ostream& operator<< (ostream& os, const Block_t& b);
		friend ostream& operator<< (ostream& os, SyTensor_t& SyT);
		friend SyTensor_t operator* (SyTensor_t& Ta, SyTensor_t& Tb);
		friend bool operator== (const Block_t& b1, const Block_t& b2);
	private:
		Qnum_t qnum;
		double* elem;
		int Rnum;		//number of rows of the block
		int Cnum;		//number of columns of the block
		int64_t offset;	//index of the first element of a block element in Tensor
};
