#include <iostream>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include "Hamiltonian.h"

int main(){
  cout<<Heisenberg_U1();
  cout<<Heisenberg(1);
  Matrix mx = matSx(1);
  cout<<mx;
  cout<<mx.eigh()[0];
  cout<<mx.eigh()[1];
}
