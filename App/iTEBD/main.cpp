#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "TensorLib.h"
#include "tools.cpp"
const int CHI = 4;

int main(){
	Qnum_t q00(0, 0);
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	for(int i = 0; i < CHI; i++)
		qnums.push_back(q00);
	
	Bond_t chi_bdr(BD_ROW, qnums);
	Bond_t chi_bdc(BD_COL, qnums);
	qnums.clear();
	qnums.push_back(q00);
	qnums.push_back(q00);
	Bond_t bdr(BD_ROW, qnums);
	Bond_t bdc(BD_COL, qnums);

	bonds.push_back(chi_bdr);
	bonds.push_back(chi_bdr);
	bonds.push_back(bdc);
	SyTensor_t Ga(bonds, "Ga");
	SyTensor_t Gb(bonds, "Gb");
	double G_elem[32];
	for(int i = 0; i < 32; i++)
		G_elem[i] = i * 0.1;
	Ga.addRawElem(G_elem);
	Gb.addRawElem(G_elem);
	cout<<Ga;
	cout<<Gb;

	bonds.clear();
	bonds.push_back(bdr);
	bonds.push_back(bdr);
	bonds.push_back(bdc);
	bonds.push_back(bdc);
	DOUBLE H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

	SyTensor_t H(bonds, "H");
	H.addRawElem(H_elem);

	bonds.clear();
	bonds.push_back(chi_bdr);
	bonds.push_back(chi_bdc);

	SyTensor_t La(bonds, "La");
	SyTensor_t Lb(bonds, "Lb");
	La.eye();
	Lb.eye();

	//cout<<H;
	//cout<<Gb;
	//cout<<La;
	SyTensor_t expH = H;
	Block_t H_blk = H.getBlock(q00);
	Block_t expH_blk = H_blk.takeExp(-1);
	expH.elemset(q00, expH_blk.getElem());
	cout<<expH;

	Network_t itebd("Itebd");
	itebd.replaceWith(0, &Ga);
	itebd.replaceWith(1, &La);
	itebd.replaceWith(2, &Gb);
	itebd.replaceWith(3, &Lb);
	itebd.replaceWith(4, &expH);
	SyTensor_t C = itebd.launch();
	//cout<<C;
	int Lb_label[] = {-1, 6};
	Lb.addLabel(Lb_label);
	SyTensor_t TH = Lb * C;
	int TH_label[] = {-1, -2, -3, -4};

	TH.reshape(TH_label, 2);
	cout << TH;
	Block_t TH_blk = TH.getBlock(q00);
	double *U = (double*) malloc(TH_blk.row() * TH_blk.row() * sizeof(DOUBLE));
	double *S = (double*) malloc(TH_blk.row() * sizeof(DOUBLE));
	double *vT = (double*) malloc(TH_blk.col() * TH_blk.col() * sizeof(DOUBLE));

	myDgesvd(TH_blk.getElem(), TH_blk.row(), TH_blk.col(), U, S, vT);
	printDMatrix(U, TH_blk.row(), TH_blk.col());
	printDMatrix(vT, TH_blk.row(), TH_blk.col());

	Block_t La_blk = La.getBlock(q00);
	double* La_elem = La_blk.getElem();
	int La_Cnum = La_blk.col();
	for(int i = 0; i < La_Cnum; i++)
		La_elem[i * La_Cnum + i] = S[i] / S[0];	//La updation over;
	
	Gb.elemset(q00, vT);
	int C_label[] = {11, 12, 13, 14};
	int Gb_label[] = {15, 13, 14};
	C.addLabel(C_label);
	Gb.addLabel(Gb_label);
	Ga = C * Gb;
	cout << Ga;
	cout << Gb;

	SyTensor_t invLb = Lb;
	Block_t invLb_blk = invLb.getBlock(q00);
	double* invLb_elem = invLb_blk.getElem();
	int invLb_Cnum = invLb_blk.col();
	for(int i = 0; i < invLb_Cnum; i++)
		invLb_elem[i * invLb_Cnum + i] = 1.0 / invLb_elem[i * invLb_Cnum + i];

	int invLb_label[] = {13, -3};
	invLb.addLabel(invLb_label);
	Gb *= invLb;	//Gb updation over
	H.check();
}
