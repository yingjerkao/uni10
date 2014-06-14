#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
void truncateU1(const UniTensor& GS, UniTensor& Al, UniTensor& Bl, Bond& bdi, Bond& bdo, Bond& bDi, Bond& bDo);
void mergeSort(vector<double>& svs, vector<Qnum>& qnums, Qnum& q, Matrix& sv_mat);
const int M = 20;

int main(){
	/*** Initialization ***/
	Qnum q0(0);
	Qnum qm1(-1);
	Qnum q1(1);
	Qnum qm2(-2);
	Qnum q2(2);
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

	int N = 20;
	int D = 2;
  Bond bDi = bdi;
  Bond bDo = bdo;
	Network HLn("HL.net");
	Network HRn("HR.net");

  qnums.clear();
  qnums.push_back(qm2);
  qnums.push_back(q0);
  qnums.push_back(q0);
  qnums.push_back(q2);
	Bond tmp_bd(BD_IN, qnums);

	vector<Bond> bond3;
	bond3.push_back(tmp_bd);
	bond3.push_back(bdo);
	bond3.push_back(bdo);
  UniTensor Al(bond3);
  UniTensor Bl(bond3);


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


		Matrix block = SB.getBlock(q0);
		vector<Matrix> rets = block.diagonalize();
		cout<<"N = "<< 2 * l<<", D = " << bDo << setprecision(10) << ", E = " << rets[0][0]  << ", e = " << rets[0][0] / (2 * l) <<endl;

		UniTensor GS(stateBonds, "GS");
		GS.elemSet(rets[1].elem());
		GS.permute(2);
		assert(GS.elemNum() == rets[1].col());
    truncateU1(GS, Al, Bl, bdi, bdo, bDi, bDo);
    //Al.printRawElem();

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

void truncateU1(const UniTensor& GS, UniTensor& Al, UniTensor& Bl, Bond& bdi, Bond& bdo, Bond& bDi, Bond& bDo){
  map<Qnum, Matrix> blocks = GS.getBlocks();
  map<int, Qnum> reorder;
  map<Qnum, vector<Matrix> > svds;
  for(map<Qnum, Matrix>::iterator it = blocks.begin(); it != blocks.end(); ++it){
    svds[it->first] = it->second.svd();
    int ord = it->first.U1() * 2;
    if(ord < 0)
      ord = (ord * -1) - 1;
    reorder[ord] = it->first;
  }
  vector<Qnum> qnums;
  vector<double> svs;
  for(map<int, Qnum>::iterator it = reorder.begin(); it != reorder.end(); ++it){
    mergeSort(svs, qnums, it->second, svds[it->second][1]);
  }

  for(int i = 0; i < qnums.size(); i++)
     cout<<qnums[i]<<", "<<setprecision(17)<<svs[i]<<endl;;
  map<Qnum, int> degs;
  std::map<Qnum,int>::iterator it;
  for(int q = 0; q < qnums.size(); q++){
    it = degs.find(qnums[q]);
    if(it != degs.end())
      it->second++;
    else
      degs[qnums[q]] = 1;
  }
  qnums.clear();
  for(it = degs.begin(); it != degs.end(); it++){
    std::map<Qnum,int>::iterator it_neg = degs.find(-it->first);
    if(it_neg == degs.end()){
      cout<<it->first<<": "<<it->second<<endl;
      assert(false);
    }
    int num = it->second > it_neg->second ? it->second : it_neg->second;
    for(int i = 0; i < num; i++)
      qnums.push_back(it->first);
  }

  //for(int i = 0; i < qnums.size(); i++)
  //   cout<<qnums[i]<<endl;


  bDi.assign(BD_IN, qnums);

  vector<Bond> bondA;
  bondA.push_back(bDi);
  bondA.push_back(bDo);
  bondA.push_back(bdo);

  vector<Bond> bondB;
  bondB.push_back(bDi);
  bondB.push_back(bdo);
  bondB.push_back(bDo);

  bDo.assign(BD_OUT, qnums);

  Al.assign(bondA);
  Bl.assign(bondB);
	vector<Qnum> blockQs = Al.blockQnum() ;
  for(int q = 0; q < blockQs.size(); q++){
    svds[blockQs[q]][0].transpose();
    Al.elemSet(blockQs[q], svds[blockQs[q]][0].elem());
    Bl.elemSet(blockQs[q], svds[blockQs[q]][2].elem());
  }
}

void mergeSort(vector<double>& svs, vector<Qnum>& qnums, Qnum& q, Matrix& sv_mat){
  if(svs.size()){
    int len = qnums.size() + sv_mat.col();
    len = len < M ? len : M;
    vector<double> ori_svs = svs;
    vector<Qnum> ori_qnums = qnums;
    Qnum q0(0);
    svs.assign(len, 0);
    qnums.assign(len, q0);
    int cnt = 0;
    int cur1 = 0;
    int cur2 = 0;
    while(cnt < len){
      if(cur1 < ori_svs.size() && cur2 < sv_mat.col()){
        if(ori_svs[cur1] >= sv_mat[cur2] - 1E-10){
          svs[cnt] = ori_svs[cur1];
          qnums[cnt] = ori_qnums[cur1];
          cur1++;
        }
        else{
          svs[cnt] = sv_mat[cur2];
          qnums[cnt] = q;
          cur2++;
        }
      }
      else if(cur2 < sv_mat.col()){
        svs[cnt] = sv_mat[cur2];
        qnums[cnt] = q;
        cur2++;
      }
      else{
        svs[cnt] = ori_svs[cur1];
        qnums[cnt] = ori_qnums[cur1];
        cur1++;
      }
      cnt++;
    }
  }
  else{
    qnums.assign(sv_mat.col(), q);
    svs.assign(sv_mat.col(), 0);
    for(int i = 0; i < sv_mat.col(); i++)
      svs[i] = sv_mat[i];
  }
}

