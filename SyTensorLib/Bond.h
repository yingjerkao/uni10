#include <iostream>
#include <iomanip>
#include <assert.h>
#include <vector>
using namespace std;
//Bond property
enum bondType{
	BD_ROW = 1,
	BD_COL = -1
};

class SyTensor_t;
class Bond_t {
	public:
		Bond_t(bondType, vector<Qnum_t>& qnums);
		Bond_t(const Bond_t& _b):type(_b.type), dim(_b.dim), Qnums(_b.Qnums), Qdegs(_b.Qdegs), offsets(_b.offsets){
			//cout<<"Copying Bond "<< this <<" from " << &_b << endl;
		}
		friend class SyTensor_t;
		friend ostream& operator<< (ostream& os, const Bond_t& b);
		friend ostream& operator<< (ostream& os, SyTensor_t& SyT);
		~Bond_t();
	private:
		bondType type;
		int dim;
		vector<Qnum_t>Qnums;	//Quantum numbers
		vector<int>Qdegs;	//Degeneracy in each quantum sector
		vector<int>offsets;	
};
