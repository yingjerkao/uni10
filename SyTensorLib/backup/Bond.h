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

class Bond_t {
	public:
		Bond_t(bondType, vector<Qnum_t>& qnums);
		friend ostream& operator<< (ostream& os, const Bond_t& b);
		~Bond_t();
	private:
		bondType type;
		int dim;
		vector<Qnum_t>Qnums;	//Quantum numbers
		vector<int>Qdegs;	//Degeneracy in each quantum sector
		vector<int>offsets;	
};
