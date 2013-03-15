#include <iostream>
#include <assert.h>
using namespace std;
class Qnum_t {
	public:
		Qnum_t()
		: U1(0), prt(0)
		{
			counter++;
			cout<<"Constructing...\n";
		}
		Qnum_t(Qnum_t& _q)
		:U1(_q.U1), prt(_q.prt)
		{
			assert(U1 < 100 && prt < 2);
			counter++;
			cout<<"Copying...\n";
		}
		~Qnum_t(){
			cout<<"DEstructing...\n";
			counter--;
		};
		const Qnum_t& operator= (const Qnum_t& _q){
			if(&_q == this)
				return *this;
			cout<<"Copying by operator=...\n";
			U1 = _q.U1, prt = _q.U1;
		}

		bool operator< (const Qnum_t& q){
			return (U1 * 2 + prt) < (q.U1 * 2 + q.prt);
		}

		int U1;
		int prt;
		static int counter;
};
int Qnum_t::counter = 0;

