#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "SyTensor.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q00(0, 0);
	Qnum_t q_10(-1, 0);
	vector<Bond_t> Abonds;		
	vector<Bond_t> Bbonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	Bond_t bd(BD_ROW, qnums);
	Abonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);
	bd.assign(BD_ROW, qnums);
	Abonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q10);
	bd.assign(BD_ROW, qnums);
	Abonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	Abonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	Abonds.push_back(bd);

	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);qnums.push_back(q10);
	bd.assign(BD_ROW, qnums);
	Bbonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	bd.assign(BD_ROW, qnums);
	Bbonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	Bbonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);
	bd.assign(BD_COL, qnums);
	Bbonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	Bbonds.push_back(bd);

	SyTensor_t SyTa(Abonds, "Ta");
	int label_tmp[] = {1, 2, 3, -1, -2};
	vector<int> labels(label_tmp, label_tmp + sizeof(label_tmp) / sizeof(int));
	SyTa.addLabel(labels);

	int tmp1[] = {-1, 0, 1, 1};
	int tmp2[] = {-1, 0, 0};
	int tmp3[] = {-1, 1};
	int tmp4[] = {-1, 0, 0, 1};
	int tmp5[] = {-1, -1, 0, 1, 1};
	int row_qnums[24];
	int cnt = 0;
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 2; k++){
				row_qnums[i * 6 + j * 2 + k] = tmp1[i] + tmp2[j] + tmp3[k];
				cnt++;
			}
	printf("cnt = %d\n", cnt);
	int col_qnums[20];
	cnt = 0;
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 5; j++){
			col_qnums[i * 5 + j] = tmp4[i] + tmp5[j];
			cnt++;
		}
	printf("cnt = %d\n", cnt);
	double elem[24 * 20];
	cnt = 1;
	for(int i = 0; i < 24; i++)
		for(int j = 0; j < 20; j++){
			if(row_qnums[i] == col_qnums[j]){
				elem[i * 20 + j] = cnt;
				cnt++;
			}
			else
				elem[i * 20 + j] = 0;
		}

	SyTa.addRawElem(elem);

	SyTensor_t SyTb(Bbonds, "Tb");
	int label2_tmp[] = {-3, -4, 1, 2, 3};
	vector<int> labels2(label2_tmp, label2_tmp + sizeof(label2_tmp) / sizeof(int));
	SyTb.addLabel(labels2);
	int label3_tmp[] = {3, -1,  2, -2, 1} ;
	vector<int> labels3(label3_tmp, label3_tmp + sizeof(label3_tmp) / sizeof(int));
	int label4_tmp[] = {2, -4, 3,  1, -3};
	vector<int> labels4(label4_tmp, label4_tmp + sizeof(label4_tmp) / sizeof(int));
	SyTa.reshape(labels3, 1);
	SyTb.reshape(labels4, 4);
			            
	cout<< SyTa;
	SyTb.randomize();
	cout<< SyTb;
	//printRawElem(SyTa);
	//printRawElem(SyTb);
	SyTensor_t SyTc = SyTa * SyTb;
	cout<<SyTc;
	SyTa.check();
}
