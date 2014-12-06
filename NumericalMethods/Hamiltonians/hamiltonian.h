#include <iostream>
#include "uni10.hpp"

void spin_check(float spin);
uni10::Matrix Sp(float spin=0.5);
uni10::Matrix Sm(float spin=0.5);
uni10::Matrix Sz(float spin=0.5);
uni10::Matrix matSx(float spin=0.5);
uni10::Bond spin_bond(float spin, uni10::bondType btype, bool U1=false);

uni10::UniTensor Heisenberg(float spin=0.5, double J=1.0);
uni10::UniTensor Heisenberg_U1(float spin=0.5, double J=1.0);

uni10::UniTensor transverseIsing(float spin, float h);
