#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;


int main(){
	Qnum q0(0);
	Qnum q1(1);
	Qnum q_1(-1);
	vector<Bond> bonds;
	vector<Bond> bonds2;
	// bond for one-in one-out tensor
	vector<Bond> two_leg;
    //
	vector<Qnum> qnums;
	qnums.push_back(q0);qnums.push_back(q0);
	Bond bdi(BD_IN, qnums);
	Bond bdo(BD_OUT, qnums);
	two_leg.push_back(bdi);
	two_leg.push_back(bdo);
	bonds.push_back(bdi);
	bonds.push_back(bdi);
	bonds.push_back(bdo);
	bonds.push_back(bdo);
	bonds2.push_back(bdi);
	bonds2.push_back(bdi);
	bonds2.push_back(bdi);
	bonds2.push_back(bdo);
	bonds2.push_back(bdo);
    bonds2.push_back(bdo);
	double H_elem[] = {1.0/4,      0,      0,     0,
						   0, -1.0/4,  1.0/2,     0,
						   0,  1.0/2, -1.0/4,     0,
						   0,      0,      0, 1.0/4};

// Spin operators 
	double S_zelement[] = { 0.5,    0,
					          0, -0.5};
	double S_pelement[] = { 0, 1.0,
		                    0,   0};
	double S_melement[] = {  0, 0,
             		       1.0, 0};

	UniTensor S_z(two_leg);
	UniTensor S_p(two_leg);
	UniTensor S_m(two_leg);
    S_z.addRawElem(S_zelement);
	S_p.addRawElem(S_pelement);
	S_m.addRawElem(S_melement);
//	dont affect them anyhow lol.
//
//Identities 
   UniTensor I_2(two_leg);
   UniTensor I_4(bonds);
   I_2.eye();
   I_4.eye();

//Labels
   int label0123[] = {0, 1, 2, 3};
   int label1032[] = {1, 0, 3, 2};
   int label4567[] = {4, 5, 6, 7};
   int label1256[] = {1, 2, 5, 6};
   int label2345[] = {2, 3, 4, 5};
   int label0_7[] = {0, 1, 2, 3, 4, 5, 6, 7};
   int label23[] = {2, 3};
   int label01[] = {0, 1};
   int label78[] = {7, 8};
// H_l - dot operators:
// ``b'' denotes boundary.
	UniTensor I_b = I_2;
	I_b.addLabel(label23);
	UniTensor S_z_b = I_b * S_z;
	UniTensor S_p_b = I_b * S_p;
	UniTensor S_m_b = I_b * S_m;

//
//
//	dot-dot operators: Hamiltonian on the left(right)-hand side: H_l(r) 

	
	UniTensor spin_z_l = S_z;
	UniTensor spin_z_r = S_z;
	//cout<< spin_z_l;
	//cout<< spin_z_r;
	UniTensor spin_p_l = S_p;
	UniTensor spin_p_r = S_p;
	UniTensor spin_m_l = S_m;
	UniTensor spin_m_r = S_m;
	//cout<< spin_z;
	int one_leg_label_1[] = {0, 2};
	int one_leg_label_2[] = {1, 3};
	
	spin_z_l.addLabel(one_leg_label_1);
	//cout<< spin_z_l;
	spin_z_r.addLabel(one_leg_label_2);
	//cout<< spin_z_r;
	spin_z_l *= spin_z_r;
	//int permute[] = {0, 1, 2, 3};
	//spin_z_l.permute(permute,2);
	//cout<< spin_z_l;
	//exit(0);
	spin_p_l.addLabel(one_leg_label_1);
	spin_m_r.addLabel(one_leg_label_2);
	spin_p_l *= spin_m_r;
	
	spin_m_l.addLabel(one_leg_label_1);
	spin_p_r.addLabel(one_leg_label_2);
	spin_m_l *= spin_p_r;
	//int two_leg_label[] = {0, 1, 2, 3};
	//spin_z_l.addLabel(two_leg_label);
	//spin_p_l.addLabel(two_leg_label);
	//spin_m_l.addLabel(two_leg_label);
	//spin_p_l.permute(permute, 2);
	//spin_m_l.permute(permute, 2);
	UniTensor dot_dot = 0.5* spin_p_l + 0.5* spin_m_l + spin_z_l;
	dot_dot.addLabel(label2345);
	//cout<< dot_dot;
	
//    
// demonstrating constructing dot-dot Hamiltonian
//---------------------------------------------------
//------------------constructing left-dot and dot-right hamiltonian-------------------------



//
//
//
//
//
//
	UniTensor H0(bonds);
	H0.addRawElem(H_elem);
	UniTensor H_dot_dot = H0;

	int l0145[] = {0, 1, 4, 5};
	H0.addLabel(l0145);
	//cout<<H0;

	UniTensor I(bonds);
	I.eye();
	int l4567[] = {2, 3, 6, 7};
	I.addLabel(l4567);
	//cout<<I;

	H0 *= I;
	//cout<<H0;
	
	int l01234567[] = {0, 1, 2, 3, 4, 5, 6, 7};
	int l20136457[] = {2, 0, 1, 3, 6, 4, 5, 7};
	int l23016745[] = {2, 3, 0, 1, 6, 7, 4, 5};
	H0.permute(l01234567, 4);
	UniTensor H1 = H0;
	H1.permute(l20136457, 4);
	UniTensor H2 = H0;
	H2.permute(l23016745, 4);
	UniTensor Hf = H0;
	Hf += H1;
	Hf += H2;
	//cout<<Hf;
	//Hf.printRawElem();
	Matrix mat = Hf.getBlock(q0);
	//cout<<mat;
   
	vector<Matrix> diag_inf;
	diag_inf = mat.diagonalize();
    //double temp[] = diag_inf[1].elem;
	Matrix state_g(1, mat.row(),diag_inf[1].elem());
	Matrix state_gd = state_g;
	state_gd.transpose();
	state_gd *= state_g;
	//cout<<"Pure Density Matrix of GS:\n";
    //cout<< state_gd;	
	//cout<< state_gd;
	UniTensor Rho = Hf;
	Rho.putBlock(q0, state_gd);
	//cout<< Rho;

	UniTensor reduced_Rho = Rho;
	reduced_Rho.partialTrace(3, 7);
	reduced_Rho.partialTrace(2, 6);
	//cout<<reduced_Rho;

	mat = reduced_Rho.getBlock(q0);
	diag_inf = mat.diagonalize();
	Matrix eigen_states = diag_inf[1];
	//cout<<eigen_states;	

	Matrix mat_HL(eigen_states.row() / 2, eigen_states.col());
	memcpy(mat_HL.elem(), (eigen_states.elem() + (eigen_states.row() - 1) * eigen_states.col()), eigen_states.col() * sizeof(double));
	memcpy(mat_HL.elem() + eigen_states.col(), (eigen_states.elem() + (eigen_states.row() - 2) * eigen_states.col()), eigen_states.col() * sizeof(double));
	//cout<<mat_HL;
	Matrix mat_HR = mat_HL;
	mat_HR.transpose();
	Matrix mat_H0(4, 4, H_elem);
	//cout<<mat_H0;
	Matrix mat_H_renormalize = mat_HL * mat_H0 * mat_HR;
	//cout<< mat_H_renormalize;
	//cout<<mat_HL * mat_H0 * mat_HR;
	UniTensor H_renormalize(two_leg);
	H_renormalize.putBlock(q0, mat_H_renormalize);
	H_renormalize.addLabel(label01);
	//cout<< H_renormalize;
	//OPerators and their renormalization

	UniTensor S_z_r(two_leg); // the subscript ''r'' denoted ''renormalized.''
	UniTensor S_p_r(two_leg);
	UniTensor S_m_r(two_leg);
	Matrix matS_z = S_z_b.getBlock(q0);
	Matrix matS_p = S_p_b.getBlock(q0);
	Matrix matS_m = S_m_b.getBlock(q0);
	//renormalization of spin operators
	Matrix matS_z_r = mat_HL * matS_z * mat_HR;
	Matrix matS_p_r = mat_HL * matS_p * mat_HR;
	Matrix matS_m_r = mat_HL * matS_m * mat_HR;
	S_z_r.putBlock(q0, matS_z_r);
	S_p_r.putBlock(q0, matS_p_r);
	S_m_r.putBlock(q0, matS_m_r);
	S_z_r.addLabel(label01);
	S_p_r.addLabel(label01);
	S_m_r.addLabel(label01);
	//^^^^ These three are renormalized operators.
	//Operators and their renormalization
	
	//boundary operators 
    //UniTensor H_boundary = H_renormalize * I_b + S_p_r	

//add labels to spin operators
	S_z.addLabel(label23);
	S_p.addLabel(label23);
	S_m.addLabel(label23);
	// include dot_dot operator
	UniTensor tmp_1 = I_2;
	tmp_1.addLabel(label78);
	UniTensor H_int = I_2 * dot_dot * tmp_1;
	H_int.addLabel(label0_7);
	//
	UniTensor H_boundary = H_renormalize * I_b + 0.5 * S_p_r * S_m + 0.5 * S_m_r * S_p + S_z_r * S_z;
	H_boundary.addLabel(label0123);
	//cout<< H_boundary;
	UniTensor H_copy = H_boundary.permute(label1032, 2);
	H_copy.addLabel(label4567);
	//cout<< H_copy;
    UniTensor H_total = H_boundary * H_copy;
	H_total.addLabel(label0_7);
	H_total += H_int;
	cout<< H_total;
	//boundary operators
	/*
	SyTensor_t H = H0;

	H.addLabel(labels1);
	I.addLabel(labels6);
	H *= I;
	I.addLabel(labels7);
	H *= I;
	H.reshape(labels_f1, 4);
	SyTensor_t Hf = H;
	H.reshape(labels_f2, 4);
	Hf += H;
	H.reshape(labels_f3, 4);
	Hf += H;
	cout<<Hf;

	*/

	//Hf.printRawElem();
	//cout<<"Trace: "<<Hf.trace()<<endl;
	//cout<<Hf.partialTrace(3, 7);

	//Hf.printRawElem();
	//cout<<"Trace: "<<Hf.trace()<<endl;
	//cout<<Hf.partialTrace(4, 8);
	//Hf.printRawElem();
	//cout<<"Trace: "<<Hf.trace()<<endl;
	//cout<<Hf;
	//SyTensor_t Hf1 = Hf;
	//cout<<Hf1;
	//cout<<Hf;
	//cout<<Hf1;
	//cout<<"Trace Hf1: "<<Hf1.trace(Hf)<<endl;
	
}

