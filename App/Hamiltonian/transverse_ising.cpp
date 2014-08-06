#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
  const double J = 1.0;
  const double g = 0.5;
  const int d = 2;
	double H_elem[] = {
     J, -g/2, -g/2,    0,
  -g/2,   -J,    0, -g/2,
  -g/2,    0,   -J, -g/2,
	 	 0, -g/2, -g/2,    J};

  // Without symmetry
	Bond bdi_d(BD_IN, d);
	Bond bdo_d(BD_OUT, d);

	vector<Bond> bond4;
	bond4.push_back(bdi_d);
	bond4.push_back(bdi_d);
	bond4.push_back(bdo_d);
	bond4.push_back(bdo_d);

  UniTensor H(bond4);
  H.addRawElem(H_elem);

  H.save("transIsing.ham");
}
