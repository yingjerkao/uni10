#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
vector<Qnum> setTruncate(vector<Bond>& bonds, map<Qnum, int>& Ms);
const int M = 20;

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Qnum qm1(-1);
	Qnum q1(1);
	Qnum qm2(-2);
	Qnum q2(2);
	Qnum qm3(-3);
	Qnum q3(3);
	Qnum qm4(-4);
	Qnum q4(4);
	Qnum qm5(-5);
	Qnum q5(5);
	vector<Qnum> qnums;
	qnums.push_back(qm1);
	qnums.push_back(q1);
	Bond bdi(BD_IN, qnums);
	Bond bdo(BD_OUT, qnums);
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

	UniTensor HL = H0;
	UniTensor HR = H0;

	int N = 30;
	int D = 2;
	Bond bDi = bdi;
	Bond bDo = bdo;
	Network HLn("HL.net");
	Network HRn("HR.net");

	map<Qnum, int> Ms;	
	Ms[q0] = 20;
	Ms[qm1] = 20;
	Ms[q1] = 20;
	Ms[qm2] = 10;
	Ms[q2] = 10;
	Ms[qm3] = 10;
	Ms[q3] = 10;
	Ms[qm4] = 10;
	Ms[q4] = 10;
	Ms[qm5] = 10;
	Ms[q5] = 10;

	/*** END initilization ***/

	for(int l = 2; l < N; l++){
		/*** Make a superblock ***/
		bond2.clear();
		bond2.push_back(bDi);
		bond2.push_back(bDo);
		ID.assign(bond2);
		ID.eye();

		UniTensor IDd(HL.bond());
		IDd.eye();
		UniTensor IdD(HR.bond());
		IdD.eye();	
		UniTensor SB = otimes(HL, IdD) + otimes(otimes(ID, H0), ID) + otimes(IDd, HR);
		//cout<<SB;
		/*** END superblock ***/
		UniTensor HL2 = otimes(HL, Id);
		HL2 += otimes(ID, H0);
		UniTensor HR2 = otimes(Id, HR);
		HR2 += otimes(H0, ID);

		vector<Bond> stateBonds;	
		stateBonds.push_back(bDi);
		stateBonds.push_back(bdi);
		stateBonds.push_back(bdi);
		stateBonds.push_back(bDi);

		vector<Bond> tmpBonds;
		tmpBonds.push_back(bDi);
		tmpBonds.push_back(bdi);
		qnums = setTruncate(tmpBonds, Ms);
		bDi.assign(BD_IN, qnums);
		vector<Bond> bond3;
		bond3.push_back(bDi);
		bond3.push_back(bDo);
		bond3.push_back(bdo);
		bDo.assign(BD_OUT, qnums);
		UniTensor Al(bond3, "Al");
		UniTensor Bl(bond3, "Bl");
		vector<Qnum> blockQs = Al.blockQnum() ;

		Matrix block = SB.getBlock(q0);
		vector<Matrix> rets = block.diagonalize();
		cout<<"N = "<< 2 * l<<", D = " << bDo << setprecision(10) << ", E = " << rets[0][0]  << ", e = " << rets[0][0] / (2 * l) <<endl;

		UniTensor GS(stateBonds, "GS");
		GS.elemSet(rets[1].elem());
		int per_label[] = {0, 1, 2, 3};
		GS.permute(per_label, 2);
		assert(GS.elemNum() == rets[1].col());
		//SB.printRawElem();
		//GS.printRawElem();

		vector<Qnum> blockQ = Al.blockQnum();
			
		/*
		Matrix GSM = GS.rawElem();//(sqrt(rets[1].col()), sqrt(rets[1].col()), rets[1].elem());
		rets = GSM.svd();
		cout<<"SVD: ";
		for(int i = 0; i < rets[1].row(); i++)
			cout<<rets[1][i]<<", ";
		cout<<endl<<endl;
		*/
		for(int q = 0; q < blockQ.size(); q++){
			rets = GS.getBlock(blockQ[q]).svd();
			cout<<"svd: ";
			for(int i = 0; i < rets[1].row(); i++)
				cout<<rets[1][i]<<", ";
			cout<<endl;
			rets[0].transpose();
			Al.elemSet(blockQ[q], rets[0].elem());
			Bl.elemSet(blockQ[q], rets[2].elem());
		}
		
		//cout<<Al;
		//cout<<Bl;
		//cout<< L;

		HLn.putTensor("Al", &Al);
		HLn.putTensor("HL2", &HL2);
		HLn.putTensorT("AlT", &Al);

		HRn.putTensor("Bl", &Bl);
		HRn.putTensor("HR2", &HR2);
		HRn.putTensorT("BlT", &Bl);

		HL = HLn.launch();
		HR = HRn.launch();
	}
}

vector<Qnum> setTruncate(vector<Bond>& bonds, map<Qnum, int>& Ms){
	vector<Qnum>qnums;
	Bond bd = combine(bonds);
	map<Qnum, int> degs = bd.degeneracy();
	for(map<Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it){
		map<Qnum, int>::const_iterator itM = Ms.find(it->first);
		if(itM != Ms.end()){
			int D = it->second < itM->second ? it->second : itM->second;
			for(int d = 0; d < D; d++)
				qnums.push_back(it->first);
		}
	}
	return qnums;
}

