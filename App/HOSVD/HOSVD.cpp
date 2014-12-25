// SIAM Journal on Matrix Analysis and Applicaions Vol. 21, No. 4, April 2000, pp. 1253-1278 p.12
#include<iostream>
#include<map>
#include<cassert>
#include<cmath>
#include<iomanip>
#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<iostream>
using namespace std;
#include"uni10.hpp"
#include<boost/filesystem.hpp>
using namespace uni10;
#define d 3
int inBondDim(UniTensor &A);
int outBondDim(UniTensor &A);
// (UniTensor A, put vector of UinBond, PATH/of/.netfile/you/want)
map<string, UniTensor> HOSVD(UniTensor &A, vector<Bond> &sbd, string PATH);
map<string, UniTensor> truncatePESS(map<string, UniTensor> &A, int CHI);
//=====================================================================================================//
UniTensor constructT(map<string, UniTensor> &hosvd, string PATH);
void writeNewTNetFile(map<string, UniTensor> &hosvd, string PATH);

int main()
{
	double M_elem[ ] = {0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740, 1.4429,
						          0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,-1.7495,
						          2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,-0.2716,
                      0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740, 1.4429,
						          0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,-1.7495,
						          2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,-0.2716,
                      0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740, 1.4429,
						          0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,-1.7495,
						          2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,-0.2716};

    const int CHI = 3;
    Matrix MA(9, 9, M_elem);
    vector<Matrix> MM = MA.svd();

    //Initialization
    Qnum q0(0);
    Bond bdi_d(BD_IN, d);
    Bond bdo_d(BD_OUT, d);
    Bond bdi_CHI(BD_IN, CHI);
    Bond bdo_CHI(BD_OUT, CHI);
    vector<Bond> Bond3;
    Bond3.push_back(bdi_d);
    Bond3.push_back(bdi_CHI);
    Bond3.push_back(bdo_d);
    Bond3.push_back(bdo_CHI);

    UniTensor A(Bond3);
    A.putBlock(q0, MA);

    map<string, UniTensor> hosvd;
    map<string, UniTensor> hosvdt;

    vector<Bond> UinBond;
    UinBond.push_back(bdi_d);
    UinBond.push_back(bdi_CHI);

    hosvd = HOSVD(A, UinBond,"core.net");

    hosvdt = truncatePESS(hosvd, CHI);

    UniTensor kk1 = constructT(hosvdt, "newT.net");
    cout << kk1 << endl;
//    for(map<string, UniTensor>::iterator it = hosvdt.begin();it != hosvdt.end(); ++it)
//        cout << it->first << endl << it->second << endl;
//    for(map<string, UniTensor>::iterator it = next(hosvd.begin(),0);it != hosvd.end(); ++it)
//        cout << it->first << endl << it->second << endl;
}

map<string, UniTensor> HOSVD(UniTensor &A, vector<Bond> &sbd, string PATH){

    //In Kagome Lattice odd labels represent virtual bonds and even labels represent physical bonds
    //Greographical Representation
    //Number of subbonds in each bond of U
    int sbn = sbd.size();
    int gi = A.bondNum() / sbn;
    assert(A.bondNum() % sbn == 0);
    //Check type of U's inBond
    for(int i = 0; i < sbd.size(); i++)
        assert(sbd[i].type() == 1);
    //clone A, A1
    UniTensor A1 = A;

    //Praparing to do svd. The fisrt step is reshaping Tensor A1 to 1st model
    A1.permute(sbn);
    int iBD = inBondDim(A1);

    // iBD is the bond dimention of inBond numbers of U
    Bond bdo_Uo(BD_OUT, iBD);

    //Initialize U's bonds
    vector<Bond> bond2(sbn+1);
    for(int i = 0; i < sbn; i++)
        bond2.at(i) = sbd[i];
    bond2.at(sbn) = bdo_Uo;

    //Tmp container include  S, U and Lambda tensors; gi Unitary tensors and Lambdas
    vector<UniTensor> hosvd(2*gi);
    Qnum q0(0);
    Matrix MA;

    //Initialize labels of A1 and core tensor
    vector<int> initLabel;
    vector<int> initLabelS;
    for(int i = 0; i < gi*sbn; i++)
        initLabel.push_back(i);
    A1.setLabel(initLabel);

    for(int i = 0; i < gi; i++)
        initLabelS.push_back(i);

    //Check whether fiting the restrict of HOSVD
    //Each bond of U, include a virturl bond and a physcial bond; Uinbond = {d, CHI};
    for(int i = 0; i < sbn; i++){
        for(int j = 0; j < gi; j++)
            assert(A1.bond()[sbn*j+i].dim() == sbd[i].dim());
    }

    //Reshap Tensor A1 to different models and do svd
    UniTensor TmpU;
    UniTensor TmpL;
    int tmp;
    vector<int> pLable;
    //Geographical
    for(int i = 0; i < gi; i++){
        //Create New Labels
        for(int j = 0; j < sbn*gi; j++){
            tmp = (sbn*i+j) % (gi*sbn);
            pLable.push_back(tmp);
        }
        A1.permute(pLable, sbn);
        MA = A1.getBlock(q0);
        TmpU.assign(bond2);
        TmpU.putBlock(q0, MA.svd()[0]);
        TmpL.assign(bond2);
        TmpL.putBlock(q0, MA.svd()[1]);
        hosvd.at(i) = TmpU;
        hosvd.at(i+gi) = TmpL;
        pLable.clear();
    }
    //recycle A1
    A1.permute(initLabel, sbn);

    //Write .netfile
    FILE* corenet = fopen(PATH.c_str(), "w");

    char TL[8];
    char TbondL[32];
    //Write Tout Labels
    for(int i = 0; i <= (gi+1); i++){
        if(i == (gi+1)){
            strcpy(TbondL, "TOUT:");
            for(int j = 0; j < gi; j++){
                int a = j + gi*sbn;
                sprintf(TL, "%d", a);
                strcat(TbondL, TL);
                if(j == 0){
                    strcat(TbondL, ";");
                }
                if(j == (gi-1))
                    fprintf(corenet, "%s\n", TbondL);
                if(j != 0  && j < (gi-1)){
                    strcat(TbondL, ",");
                }
            }
        }
        if(i == gi){
            strcpy(TbondL, "A:");
            for(int j = 0; j < gi*sbn; j++){
                sprintf(TL, "%d", j);
                strcat(TbondL, TL);
                if(j == (A1.inBondNum()-1))
                    strcat(TbondL, ";");
                if(j == (gi*sbn-1))
                    fprintf(corenet, "%s\n", TbondL);
                if(j != (A1.inBondNum()-1)  && j < (gi*sbn-1))
                    strcat(TbondL, ",");
            }
        }
        if(i < gi){
            sprintf(TbondL, "U%d:", i);
            for(int j = 0; j <= sbn; j++){
                int a = i*sbn+j;
                int b = i+sbn*gi;
                if(j == sbn-1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL,";");
                }
                if(j < sbn -1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL, ",");
                }
                if(j == sbn){
                    sprintf(TL, "%d", b);
                    strcat(TbondL, TL);
                    fprintf(corenet,"%s\n", TbondL);
                }
            }
        }
    }
    fclose(corenet);

    //Get coretensor
    Network CoreTensor(PATH);
    char TUname[8];
    char TSname[8];
    char LDname[8];
    for(int i = 0; i <= gi; i++){
        if(i == gi)
            CoreTensor.putTensor("A", &A1);
        if(i < gi){
            sprintf(TUname, "U%d", i);
            CoreTensor.putTensor(TUname, &hosvd[i]);
            }
    }
    //core tensor of mode 0 --> S0
    UniTensor core = CoreTensor.launch();
    core.setLabel(initLabelS);
/**************** Reshape core tensor S ****************/
/*  vector<int> SLabel;
    for(int i = 0; i < gi; i++){
        for(int j = 0; j < gi*sbn; j++){
           tmp = (i*sbn+j)%(gi*sbn);
           SLabel.push_back(tmp);
        }
        core.permute(SLabel, sbn);
        hosvd.at(i+gi) = core;
        SLabel.clear();
    }                                                   */
/********************************************************/
    map<string, UniTensor> hosvdm;
    hosvdm["S0"] = core;
    for(int i = 0; i < gi; i++){
        sprintf(TUname, "U%d", i);
        sprintf(LDname, "L%d", i);
        hosvdm[TUname] = hosvd[i];
        hosvdm[LDname] = hosvd[i+gi];
    }
    return hosvdm;
}

map<string, UniTensor> truncatePESS(map<string, UniTensor> &A, int CHI){
    assert(A.size() % 2 == 1);
    assert(CHI <= inBondDim(A["U0"]));
    for(int i = 0; i < (A.size()-1)/2; i++){
        char TUname[8];
        sprintf(TUname, "U%d", i);
        assert( outBondDim(A[TUname]) == inBondDim(A["S0"]));
    }
    //clone A["S0"] for keeping the values
    //clone A["S)"]
    UniTensor TmpSS = A["S0"];
    map<string, UniTensor> hosvdt;

    Qnum q0(0);
    int Sbn = A["S0"].bondNum(); //bond number of core tensor
    int Ubn = A["U0"].bondNum();
    vector<int> SLabel, S1Label;
    vector<int> initLabel = A["S0"].label();

    //bond with dimension chi
    Bond bdi_x(BD_IN, CHI);
    Bond bdo_x(BD_OUT, CHI);

    vector<Bond> bondS;
    vector<Bond> bondU = A["U0"].bond();
    UniTensor TmpS;
    UniTensor TmpU;

    char Tname[8];
    int Tmp;
    Matrix MU(CHI,inBondDim(A["U0"]));
    Matrix TmpM;

    for(int i = 0; i < Sbn; i++){

        //Truncate U
        sprintf(Tname, "U%d", i);
        TmpM = MU;
        TmpM.setElem(A[Tname].getBlock(q0).transpose().getElem());
        bondU.at(Ubn-1) = bdo_x;
        TmpU.assign(bondU);
        TmpU.putBlock(q0, TmpM.transpose());
        hosvdt[Tname] = TmpU;

        //Truncate S
        //Reshaping to different modes
        for(int j = 0; j < Sbn; j++){
            Tmp = (i+j) % Sbn;
            SLabel.push_back(Tmp);
        }
        A["S0"].permute(SLabel, 1);
        TmpS = A["S0"];
        S1Label = A["S0"].label();
        bondS = A["S0"].bond();
        bondS.at(0) = bdi_x;
        A["S0"].assign(bondS);
        Matrix MS(CHI, outBondDim(TmpS));
        TmpM = MS;
        TmpM.setElem(TmpS.getBlock(q0).getElem());
        A["S0"].putBlock(q0, TmpM);
        A["S0"].setLabel(S1Label);
        SLabel.clear();
    }
    hosvdt["S0"] = A["S0"].permute(initLabel, 1);
    A["S0"] = TmpSS;
    return hosvdt;
}

int inBondDim(UniTensor &A){

    int iBD = 1;
    for(int i = 0; i < A.inBondNum(); i++)
        iBD *= A.bond()[i].dim();

    return iBD;
}

int outBondDim(UniTensor &A){

    int iBD = 1;
    for(int i = 0; i < A.bondNum()-A.inBondNum(); i++)
        iBD *= A.bond()[A.inBondNum()+i].dim();

    return iBD;
}
//=======================================================================================//
UniTensor constructT(map<string, UniTensor> &hosvd, string PATH){

  //Write .netfile
  writeNewTNetFile(hosvd, PATH);

  int gi = hosvd["S0"].bondNum();

  Network contractNewT(PATH);

  for(map<string, UniTensor>::iterator it = hosvd.find("S0"); it != hosvd.end(); ++it)
    contractNewT.putTensor(it->first, &it->second);

  UniTensor newT = contractNewT.launch();

  return newT;
}

void writeNewTNetFile(map<string, UniTensor> &hosvd, string PATH){

    int gr = hosvd["S0"].bondNum();
    int sbn = hosvd["U0"].bondNum()-1;

    FILE* newTnet = fopen(PATH.c_str(), "w");

    char TL[8];
    char TbondL[32];
    //Write Tout Labels
    for(int i = 0; i <= (gr+1); i++){
        if(i == gr){
            strcpy(TbondL, "S0:");
            for(int j = 0; j < gr; j++){
                int a = j + gr*sbn;
                sprintf(TL, "%d", a);
                strcat(TbondL, TL);
                if(j == 0){
                    strcat(TbondL, ";");
                }
                if(j == (gr-1))
                    fprintf(newTnet, "%s\n", TbondL);
                if(j != 0  && j < (gr-1)){
                    strcat(TbondL, ",");
                }
            }
        }
        if(i == (gr+1)){
            strcpy(TbondL, "TOUT:");
            for(int j = 0; j < gr*sbn; j++){
                sprintf(TL, "%d", j);
                strcat(TbondL, TL);
                if(j == (sbn-1))
                    strcat(TbondL, ";");
                if(j == (gr*sbn-1))
                    fprintf(newTnet, "%s\n", TbondL);
                if(j < (gr*sbn-1) && j !=(sbn-1))
                    strcat(TbondL, ",");
            }
        }
        if(i < gr){
            sprintf(TbondL, "U%d:", i);
            for(int j = 0; j <= sbn; j++){
                int a = i*sbn+j;
                int b = i+sbn*gr;
                if(j == sbn-1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL,";");
                }
                if(j < sbn -1){
                    sprintf(TL, "%d", a);
                    strcat(TbondL, TL);
                    strcat(TbondL, ",");
                }
                if(j == sbn){
                    sprintf(TL, "%d", b);
                    strcat(TbondL, TL);
                    fprintf(newTnet,"%s\n", TbondL);
                }
            }
        }
    }
    fclose(newTnet);
}
/* 	double M_elem[ ] = { 0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
					      	     0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           0.8924,-0.4898, 2.4288, 1.7753,-1.5077, 4.0337,-0.6631, 1.9103,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335,
	                     0.9073, 0.7158,-0.3698, 1.7842, 1.6970, 0.0151, 2.1236,-0.0740,
						           2.1488, 0.3054, 2.3753, 4.2495, 0.3207, 4.7146, 1.8260, 2.1335};*/
