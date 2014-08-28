#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
#include <time.h>
using namespace uni10;

const int M = 20;
UniTensor mirror(const UniTensor& T);

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
	Id.identity();
	UniTensor ID(bond2);
	ID.identity();

	UniTensor HL = H0;
	UniTensor HR = mirror(H0);

	int N = 20;
	int D = 2;
	Bond bDi = bdi;
	Bond bDo = bdo;
	Network HLn("HL.net");
	Network HRn("HR.net");

  Matrix psi;
  /*** END initilization ***/

	clock_t total_t = 0;
	for(int l = 2; l < N; l++){
		/*** Make a superblock ***/
		bond2.clear();
		bond2.push_back(bDi);
		bond2.push_back(bDo);
		ID.assign(bond2);
		ID.identity();

		UniTensor IDd(HL.bond());
		IDd.identity();
		UniTensor IdD(HR.bond());
		IdD.identity();
		UniTensor SB = otimes(HL, IdD) + otimes(otimes(ID, H0), ID) + otimes(IDd, HR);
		/*** END superblock ***/
		UniTensor HL2 = otimes(HL, Id);
		HL2 += otimes(ID, H0);
		UniTensor HR2 = otimes(Id, HR);
		HR2 += otimes(H0, ID);

    /*
		vector<Matrix> rets = SB.getBlock(q0).diagonalize();
		cout<<"N = "<< 2 * l << ", D = " << D << setprecision(10) << ", E = " << rets[0][0]  << ", e = " << rets[0][0] / (2 * l) <<endl;
		Matrix GS(sqrt(rets[1].col()), sqrt(rets[1].col()), rets[1].getElem());
    */

    vector<Matrix> rets;
    Matrix blk = SB.getBlock(q0);
    psi.randomize();
    int iter = 200;
	  clock_t t;
  	t = clock();
    if(psi.col() != blk.col()){
      psi.resize(1, blk.col());
      psi.randomize();
    }

    double E0 = blk.lanczosEig(psi, iter);
	  total_t += clock() - t;
		cout<<"N = "<< 2 * l<<", D = " << D << setprecision(10) << ", E = " << E0  << ", e = " << E0 / (2 * l) <<", iter = "<<iter<<endl;
		Matrix GS(sqrt(psi.col()), sqrt(psi.col()), psi.getElem());

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
		Al.addRawElem(rets[0].getElem());
		UniTensor Bl(bondB, "Bl");
		Bl.addRawElem(rets[2].getElem());

		HLn.putTensor("Al", &Al);
		HLn.putTensor("HL2", &HL2);
		HLn.putTensorT("AlT", &Al);

		HRn.putTensor("Bl", &Bl);
		HRn.putTensor("HR2", &HR2);
		HRn.putTensorT("BlT", &Bl);

		HL = HLn.launch();
		HR = HRn.launch();
	}
	printf ("It took %lu clicks (%f seconds).\n",total_t,((float)total_t)/CLOCKS_PER_SEC);
}

UniTensor mirror(const UniTensor& T){
	vector<int> labels = T.label();
	vector<int> mlabels(labels.size());
	for(int l = 0; l < labels.size(); l++){
		if(l < T.inBondNum())
			mlabels[l] = labels[T.inBondNum() - l - 1];
		else
			mlabels[l] = labels[T.inBondNum() + T.bondNum() - l - 1];
	}
	UniTensor mT = T;
	mT.permute(mlabels, T.inBondNum());
	return mT;
}
