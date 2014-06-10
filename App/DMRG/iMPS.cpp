#include <iostream>
#include <sstream>
#include <string>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;

const int M = 20;
const int N = 10;

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Bond bdi(BD_IN, 2);
	Bond bdo(BD_OUT, 2);
	vector<Bond> bond4;
	bond4.push_back(bdi);
	bond4.push_back(bdi);
	bond4.push_back(bdo);
	bond4.push_back(bdo);
	vector<Bond> bond2;
	bond2.push_back(bdi);
	bond2.push_back(bdo);


	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	UniTensor H0(bond4, "H0");
	H0.addRawElem(H_elem);

	UniTensor Id(bond2, "Id");
	Id.eye();
	UniTensor ID(bond2);
	ID.eye();

  vector<UniTensor> HLs;
  vector<UniTensor> HRs;
	HLs.push_back(H0);
	HRs.push_back(H0);
  vector<UniTensor> As;
  vector<UniTensor> Bs;
  As.push_back(Id);
  Bs.push_back(Id);
	//cout<<HL;
	//cout<<HR;

	int D = 2;
	Bond bDi = bdi;
	Bond bDo = bdo;
	Network HLn("HL.net");
	Network HRn("HR.net");

	/*** END initilization ***/

	for(int l = 1; l < N; l++){
		/*** Make a superblock ***/
		bond2.clear();
		bond2.push_back(bDi);
		bond2.push_back(bDo);
		ID.assign(bond2);
		ID.eye();

		UniTensor IDd(HLs[l - 1].bond());
		IDd.eye();
		UniTensor IdD(HRs[l - 1].bond());
		IdD.eye();
		UniTensor SB = otimes(HLs[l - 1], IdD) + otimes(otimes(ID, H0), ID) + otimes(IDd, HRs[l - 1]);
		/*** END superblock ***/
		UniTensor HL2 = otimes(HLs[l - 1], Id);
		HL2 += otimes(ID, H0);
		UniTensor HR2 = otimes(Id, HRs[l - 1]);
		HR2 += otimes(H0, ID);

		vector<Matrix> rets = SB.getBlock(q0).diagonalize();
		cout<<"N = "<< 2 * (l + 1)<<", D = " << D << setprecision(10) << ", E = " << rets[0][0]  << ", e = " << rets[0][0] / (2 * (l + 1)) <<endl;
		Matrix GS(sqrt(rets[1].col()), sqrt(rets[1].col()), rets[1].elem());

		D = GS.row();
		D = D < M ? D : M;
		bDi.assign(BD_IN, D);

		vector<Bond> bondA;
		bondA.push_back(bDi);
		bondA.push_back(bDo);
		bondA.push_back(bdo);
		vector<Bond> bondB;
		bondB.push_back(bDi);
		bondB.push_back(bdo);
		bondB.push_back(bDo);
		bDo.assign(BD_OUT, D);

		rets = GS.svd();

		UniTensor Al(bondA, "Al");
		rets[0].transpose();
		Al.addRawElem(rets[0].elem());
		UniTensor Bl(bondB, "Bl");
		Bl.addRawElem(rets[2].elem());
    As.push_back(Al);
    Bs.push_back(Bl);

		HLn.putTensor("Al", &Al);
		HLn.putTensor("HL2", &HL2);
		HLn.putTensorT("AlT", &Al);

		HRn.putTensor("Bl", &Bl);
		HRn.putTensor("HR2", &HR2);
		HRn.putTensorT("BlT", &Bl);

    HLs.push_back(HLn.launch());
    HRs.push_back(HRn.launch());
	}
  for(int l = 0; l < As.size(); l++){
    ostringstream Al;
    ostringstream Bl;
    ostringstream HRl;
    ostringstream HLl;
    string dir("Data/");
    Al<<dir<<"A"<<l;
    Bl<<dir<<"B"<<l;
    HLl<<dir<<"HL"<<l;
    HRl<<dir<<"HR"<<l;
    As[l].save(Al.str());
    Bs[l].save(Bl.str());
    HLs[l].save(HLl.str());
    HRs[l].save(HRl.str());
  }
}

