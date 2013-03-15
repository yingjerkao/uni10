#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <string.h>
#include <assert.h>
#include <stdint.h>
using namespace std;
size_t MEM;
const int INIT = 1;		//initialized
const int HAVELABEL = 2;		//tensor with label
const int HAVEELEM = 4;		//tensor with elements
const int DISPOSABLE = 8;		//The elements of tensor is disposable through 'clone', 'reshapeClone', 'reshapeElem'. The elements of disposable Tensor is read-only.
const int ELEMFREED = 16;		//the memory space of elements is freed
#define DOUBLE		double

enum bondType{
	BD_ROW = 1,
	BD_COL = -1
};
class Bond_t{
	public:
		Bond_t()
		: U1(0), prt(0)
		{
			cout<<"Constructing...\n";
		}
		Bond_t(Bond_t& _q)
		:U1(_q.U1), prt(_q.prt)
		{
			cout<<"Copying...\n";
		}
		~Bond_t(){
			cout<<"DEstructing...\n";
		};
		int U1;
		int prt;
}

