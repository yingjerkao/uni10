using namespace std;
#include <vector>
typedef struct{
	int tenNum;			//number of tensor
	vector<int> tenBond;
	vector<int*> laList;	//the label array of each tensor.
	int* outL;	//Tout label.
	vector<int> order;		//order of contraction when constructing a operator
}Diagram;


vector<Diagram*> Diagram_List;	//Read only diagram list for a entire MERA
//Setting the Diagram_List of from file.
void initNetwork(char* DiagramFile);

//Give a diagram and a series of tensors, doing the series contraction.

void operate(int DiagramId, vector<Tensor*>& Tlist, Tensor* Tout);


void initNetwork(char* fn){
	FILE *fp = fopen(fn, "r");
	assert(fp != NULL);
	Diagram* D;
	char tmp[32];
	for(int d = 0; d < DIAG_NUM; d++){
		fscanf(fp, "%s", tmp);	//read "Did:"
		fscanf(fp, "%s", tmp);	//read Did
		D = (Diagram*)calloc(1, sizeof(Diagram));
		fscanf(fp, "%d", &(D->tenNum));
		int i, j;
		int bond;
		for(i = 0; i < D->tenNum; i++){
			fscanf(fp, "%d", &bond);
			D->tenBond.push_back(bond);
		}
		int outBonds;
		fscanf(fp, "%d", &outBonds);
		int *temp;
		for(i = 0; i < D->tenNum; i++){
			temp = (int*)calloc(D->tenBond[i], sizeof(int));
			for(j = 0; j<D->tenBond[i]; j++)
				fscanf(fp, "%d", &(temp[j]));
			D->laList.push_back(temp);
		}
		vector<int> tmp_outL;
		int ol;
		for(i = 0; i < outBonds; i++){
			fscanf(fp, "%d", &ol);
			tmp_outL.push_back(ol);
		}
		D->outL = (int*)malloc(tmp_outL.size() * sizeof(int));
		for(i = 0; i < tmp_outL.size(); i ++)
			D->outL[i] = tmp_outL[i];
		int order;
		for(i = 0; i<(D->tenNum-1)*2; i++){
			fscanf(fp, "%d", &order);
			D->order.push_back(order);
		}
		Diagram_List.push_back(D);
	}
	fclose(fp);
}

void operate(int Did, vector<Tensor*>& Tlist, Tensor* Tout){
	Diagram* D = Diagram_List[Did];
	assert(D->tenNum == Tlist.size());
	int i, j; 
	// Add labels of tensors in AscendL
	for(i = 0; i < D->tenNum; i++)
		addLabel(Tlist[i], D->laList[i]);
	//Do tenNum-1 times tensor contraction, since there are tenNum's tesnsors to be contracted to one tensor
	Tensor *Tc;
	for(i = 0; i < 2*(D->tenNum-1) - 2; i += 2){
		Tc = (Tensor*)calloc(1, sizeof(Tensor));
		contraction(Tlist[D->order[i]], Tlist[D->order[i+1]], Tc);
		Tc->status |= DISPOSABLE;
		for(j = i; j < i + 2; j++)
			if((Tlist[D->order[j]]->status) & DISPOSABLE)
				recycle(Tlist[D->order[j]]);
		Tlist.erase(Tlist.begin()+D->order[i]);
		Tlist.erase(Tlist.begin()+D->order[i+1]-1);	//-1 since there are already one element before order[i+1] being erased
		Tlist.push_back(Tc);
	}
	contraction(Tlist[D->order[i]], Tlist[D->order[i+1]], Tout, D->outL);
	for(j = i; j < i + 2; j++)
		if((Tlist[D->order[j]]->status) & DISPOSABLE)
			recycle(Tlist[D->order[j]]);
}
