#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "uni10.hpp"
using namespace uni10;
#include <time.h>
#include <vector>

int main(){
  vector<Bond> bonds;
  bonds.push_back(Bond(BD_IN, 5));
  bonds.push_back(Bond(BD_OUT, 3));
  bonds.push_back(Bond(BD_OUT, 4));

  vector<Qnum> qnums;
  qnums.push_back(Qnum(-1));
  qnums.push_back(Qnum(0));
  qnums.push_back(Qnum(0));
  qnums.push_back(Qnum(1));
  vector<Qnum> qnum1s(qnums);
  qnum1s.push_back(Qnum(2));
  qnum1s.push_back(Qnum(2));
  qnum1s.push_back(Qnum(1));
  qnum1s.push_back(Qnum(1));
  qnum1s.push_back(Qnum(-1));
  qnum1s.push_back(Qnum(-1));
  vector<Bond> bsys;
  bsys.push_back(Bond(BD_IN, qnum1s));
  qnums.push_back(Qnum(-1));
  bsys.push_back(Bond(BD_OUT, qnums));
  qnums.push_back(Qnum(1));
  bsys.push_back(Bond(BD_OUT, qnums));

  cout<<"============ CUniTensor::CuniTensor() ============\n";
  CUniTensor T0;
  //cout<<T0;
  cout<<"Okay!\n";

  cout<<"============ CUniTensor(std::complex<double> val) ============\n";
  CUniTensor T1(complex<double>(3.5, 2.1));
  //cout<<T1;
  cout<<"Okay!\n";

  cout<<"============ CUniTensor(const std::vector<Bond>& _bonds, const std::string& _name = "") ============\n";
  CUniTensor T2(bonds, "T2");
  cout<<T2;
  cout<<"Okay!\n";

  cout<<"============ CUniTensor(const std::vector<Bond>& _bonds, int* labels, const std::string& _name = "") ============\n";
  int label[] = {-3, -4, -5};
  CUniTensor T3(bsys, label, "T3");
  cout<<T3;
  cout<<"Okay!\n";

  cout<<"============ randomize() ============\n";
  T2.randomize();
  T3.randomize();
  cout<<T2<<T3;
  cout<<"Okay!\n";

  cout<<"============ CUniTensor(const CUniTensor& UniT) ============\n";
  CUniTensor T4(T2);
  cout<<T4;
  cout<<"Okay!\n";

  cout<<"============ operator=(const CUniTensor& UniT) ============\n";
  T4 = T3;
  cout<<T4;
  cout<<"Okay!\n";

  cout<<"============ CUniTensor(const CBlock& UniT) ============\n";
  CMatrix M5(3, 4);
  M5.randomize();
  CUniTensor T5 = M5;
  cout<<T5;
  cout<<"Okay!\n";

  cout<<"============ assign(const std::vector<Bond>& _bond) ============\n";
  T4.assign(bonds);
  cout<<T4;
  cout<<"Okay!\n";

  cout<<"============ getRawElem(), setRawElem(), printRawElem() ============\n";
  CMatrix M3 = T3.getRawElem();
  cout<<M3;
  M3 *= 10;
  T3.setRawElem(M3);
  cout<<T3;
  cout<<"Okay!\n";

  cout<<"============ getBlock(), putBlock(), const_getBlock() ============\n";
  cout<<T3.const_getBlock(qnums[0]);
  M3 = T3.getBlock();
  M3 *= 0.1;
  T3.putBlock(M3);
  cout<<T3;
  cout<<"Okay!\n";

  cout<<"============ permute(), transpose() ============\n";
  CUniTensor tmp_T3 = T3;
  int per1[] = {-4, -5, -3};
  T3.permute(per1, 3);
  T3.permute(per1, 0);
  T3.permute(label, 1);
  cout<<(T3.elemCmp(tmp_T3))<<endl;
  T3.transpose();
  cout<<T3;
  T3.cTranspose();
  cout<<T3;
  cout<<"Okay!\n";

  cout<<"=========== Conversion ===========\n";
  UniTensor nT1(bsys, label);
  nT1.randomize();
  CUniTensor T6(nT1);
  cout<<nT1<<T6;
  T6.randomize();
  UniTensor nT2 = T6;
  cout<<nT2<<T6;

  Matrix M1(5, 3);
  Matrix M2(5, 6, true);
  CMatrix M6(4, 3);
  M1.randomize();
  M2.randomize();
  M6.randomize();
  CUniTensor T7 = M1; //Matrix to CUniTensor
  CUniTensor T8 = M2; //Diagonal matrix to CUniTensor
  UniTensor nT3(M6);  //CMatrix to UniTensor
  cout<<M1<<T7;
  cout<<M2<<T8;
  cout<<M6<<nT3;
  Matrix rm1 = nT1.getRawElem();
  CMatrix rm2 = T6.getRawElem();
  cout<<nT1<<T6;
  nT1.setRawElem(rm2); //UniTensor::setRawElem(CBlock&)
  T6.setRawElem(rm1); //UniTensor::setRawElem(CBlock&)
  cout<<nT1<<T6;

  cout<<T6 * complex<double>(3, 4);
  cout<<complex<double>(3, 4) * T6;
  cout<<nT1 * complex<double>(3, 4);
  cout<<complex<double>(3, 4) * nT1;

  cout<<nT1 + T6;
  cout<<T6 + nT1;
  cout<<nT1 * T6;
  cout<<T6 * nT1;
  cout<<(T6+=nT1);
  cout<<(nT1+=T6);


  cout<<"Okay!\n";

}
