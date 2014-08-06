#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
  const int d = 2;
	double H_elem[] = {
    1.0/4,      0,      0,     0,
        0, -1.0/4,  1.0/2,     0,
        0,  1.0/2, -1.0/4,     0,
        0,      0,      0, 1.0/4};

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

  H.save("heisenberg.ham");
  cout<<H;

  Qnum q1(1);
  Qnum qm1(-1);
  vector<Qnum> qnums(d, q1);
  qnums[1] = qm1;
  Bond bdi_dU1(BD_IN, qnums);
  Bond bdo_dU1(BD_OUT, qnums);
  bond4[0] = bdi_dU1;
  bond4[1] = bdi_dU1;
  bond4[2] = bdo_dU1;
  bond4[3] = bdo_dU1;
  UniTensor H_U1(bond4);
  H_U1.addRawElem(H_elem);
  H_U1.save("heisenberg_U1.ham");
  cout<<H_U1;
  H_U1.printRawElem();
}
