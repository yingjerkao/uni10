#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>

int main(){
	double H3_elem[] = {\
     1, 0, 0, 0, 0, 0, 0, 0, 0,\
		 0, 0, 0, 1, 0, 0, 0, 0, 0,\
		 0, 0,-1, 0, 1, 0, 0, 0, 0,\
		 0, 1, 0, 0, 0, 0, 0, 0, 0,\
		 0, 0, 1, 0, 0, 0, 1, 0, 0,\
		 0, 0, 0, 0, 0, 0, 0, 1, 0,\
		 0, 0, 0, 0, 1, 0,-1, 0, 0,\
		 0, 0, 0, 0, 0, 1, 0, 0, 0,\
		 0, 0, 0, 0, 0, 0, 0, 0, 1\
		};
  Matrix H3(9, 9, H3_elem);
  //cout<<H3;
  vector<Matrix>rets = H3.diagonalize();
  //cout<<rets[0];
  //cout<<rets[1];
  /*
  double H2_elem[] = {\
  	0.25, 0, 0, 0,\
	  0, -0.25, 0.5, 0,\
  	0, 0.5, -0.25, 0,\
  	0, 0, 0, 0.25\
  };
  Matrix H2(4, 4, H2_elem);
  rets = H2.diagonalize();
  cout<<rets[0];
  cout<<rets[1];
  */
  Matrix psi(1, 9);
  double psi_elem[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1};
  psi.setElem(psi_elem);
  int iter = 4;
  cout<<H3.lanczosEig(psi, iter, 5E-15);
  cout<<psi;
  cout<<iter<<endl;


}

