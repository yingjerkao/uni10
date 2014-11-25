#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <time.h>
using namespace std;
#include "uni10.hpp"
using namespace uni10;

#include "iTEBD_U1_tools.cpp"
#include "hamiltonian.cpp"


int main(){
  // Define the parameters of the model / simulation
  int chi = 40;
  double delta = 0.1;
  int N = 2000;
  UniTensor H = Heisenberg_U1();

  Bond bdi_mid = H.bond(0);
  Bond bdo_mid = H.bond(2);
  bdi_mid.combine(H.bond(0));
  bdo_mid.combine(H.bond(2));

  vector<UniTensor> Gs;
	vector<Bond> bondG;
	bondG.push_back(H.bond(0));
	bondG.push_back(bdo_mid);
	bondG.push_back(H.bond(2));
  Gs.push_back(UniTensor(bondG, "GA"));

  bondG[0] = bdi_mid;
  bondG[1] = H.bond(2);
  Gs.push_back(UniTensor(bondG, "GB"));
  Gs[0].randomize(), Gs[1].randomize();

  vector<UniTensor> Ls;
  vector<Bond> bondL;
  bondL.push_back(bdi_mid);
  bondL.push_back(bdo_mid);
  Ls.push_back(UniTensor(bondL, "LA"));

  bondL[0] = H.bond(0);
  bondL[1] = H.bond(2);
  Ls.push_back(UniTensor(bondL, "LB"));
  Ls[0].randomize(), Ls[1].randomize();

  UniTensor U(H.bond());
  vector<Qnum> blk_qnums = H.blockQnum();
  for(vector<Qnum>::iterator it = blk_qnums.begin(); it != blk_qnums.end(); it++)
    U.putBlock(*it, takeExp(-delta, H.getBlock(*it)));

  int ordA[] = {-1, 3, 1};
  int ordB[] = {3, -3, 2};
  int ordU[] = {1, 2, -2, -4};

  UniTensor theta;
  // Perform the imaginary time evolution alternating on A and B bonds
  for(int step = 0; step < N; step++){
    // Construct theta
    if(step % 100 == 0)
      cout<<" ========================== Step "<<step <<" =======================\n";
    int A = step % 2;
    int B = (step + 1) % 2;
    bondcat(Gs[A], Ls[A], 1);
    bondcat(Gs[A], Ls[B], 0);
    bondcat(Gs[B], Ls[B], 1);
    Gs[A].setLabel(ordA);
		Gs[B].setLabel(ordB);
		U.setLabel(ordU);
		theta = contract(Gs[A], Gs[B], true);
		theta *= U;
		int ordTheta[] = {-1, -2, -3, -4};
		theta.permute(ordTheta, 2);
    // update Gs[A], Gs[B], Ls[A]
    setTruncation(theta, Gs[A], Gs[B], Ls[A], chi);
		Gs[A].permute(ordA, 1);
    bondrm(Gs[A], Ls[B], 0);
    bondrm(Gs[B], Ls[B], 1);
  }
  UniTensor val = theta * theta;
  cout<<"E = "<<setprecision(12)<<-log(val[0])/delta/2<<endl;
}
