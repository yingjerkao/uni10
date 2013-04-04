#include <iostream>
#include <assert.h>
#include <map>
using namespace std;
#include "SyTensor.h"

int main(){
	Qnum_t q10(1, 0);
	Qnum_t q00(0, 0);
	Qnum_t q_10(-1, 0);
	vector<Bond_t> bonds;		
	vector<Qnum_t> qnums;
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	Bond_t bd(BD_ROW, qnums);
	bonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);
	bd.assign(BD_ROW, qnums);
	bonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q10);
	bd.assign(BD_ROW, qnums);
	bonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q00);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	bonds.push_back(bd);
	qnums.clear();
	qnums.push_back(q_10);qnums.push_back(q_10);qnums.push_back(q00);qnums.push_back(q10);qnums.push_back(q10);
	bd.assign(BD_COL, qnums);
	bonds.push_back(bd);

	SyTensor_t SyT(bonds, "HelloD");
	int label_tmp[] = {-1, -2, -3, -4, -5};
	vector<int> labels(label_tmp, label_tmp + sizeof(label_tmp) / sizeof(int));
	SyT.addLabel(labels);

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

	SyT.addRawElem(elem);

	SyTensor_t rSyT = SyT;
	//int label2_tmp[] = {-1, -2, -3, -4, -5};
	int label2_tmp[] = {-4, -5, -1, -2, -3};
	//int label2_tmp[] = {-2, -5, -3, -1, -4};
	vector<int> labels2(label2_tmp, label2_tmp + sizeof(label2_tmp) / sizeof(int));
	//rSyT.reshape(labels2, 2);
	rSyT.transpose();
	rSyT.eye();
			            
	cout<< SyT;
	cout<< rSyT;
	printRawElem(SyT);
	//printRawElem(rSyT);
	SyT.check();
	Block_t blk = SyT.getBlock(q10);
	cout<<blk.row()<<endl;
	cout<<blk.col()<<endl;
}

