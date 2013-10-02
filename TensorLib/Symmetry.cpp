typedef struct{
	int row_begin;
	int col_begin;
	int row;
	int col;
}Group;

typedef struct{
	vector<Group> group;
	int totRow;
	int totCol;
	//"table" is used to do basis transformation from original to block-diagonal basis
	//in our case, shaft the position of column and row, 
	int* table;
	int* invTable;
	double* factor;
	double* invFactor;
	int mapNum;
}Symmetry;


void addGroup(Symmetry* S, int row, int col, int pos = -1){
	assert(S->mapNum == 0);	//if add Transformation matrix, you CANNOT add more group.
	Group g;
	g.row = row;
	g.col = col;
	if(S->group.size()){
		g.row_begin = S->totRow;
		g.col_begin = S->totCol;
		S->totRow += row;
		S->totCol += col;
	}
	else{	//First add, do initialization
		g.row_begin = 0;
		g.col_begin = 0;
		S->totRow = row;
		S->totCol = col;
		S->table = NULL;
		S->factor = NULL;
		S->invTable = NULL;
		S->invFactor = NULL;
		S->mapNum = 0;
	}
	vector<Group>::iterator it = S->group.begin();
	if(pos == -1)
		S->group.push_back(g);
	else if(pos < S->group.size())
		S->group.insert(it + pos, g);
}
/*void addTable(Symmetry* S, int mapNum, int* table, double* factor){
	assert(S->group.size());
	int tableSz = S->totRow > S->totCol ? S->totRow: S->totCol;	//max(S->totRow, S->totCol)
	S->mapNum = mapNum;
	if(S->table != NULL)
		free(S->table);
	S->table = (int*)calloc(mapNum * tableSz, sizeof(int));
	memcpy(S->table, table, sizeof(int)*(mapNum * tableSz));
	if(S->invTable != NULL)
		free(S->invTable);
	S->invTable = (int*)calloc(mapNum * tableSz, sizeof(int));
    for (int i=0; i< mapNum * tableSz; i++)
    	S->invTable[i]  = 0;
	
    #ifdef SANITY_CHECK
	int check[tableSz];
	memset(check, 0, tableSz * sizeof(int));
	for(int i = 0; i < tableSz * mapNum; i++){
		if(mapNum == 1)
			check[table[i]] ++;
		else if(factor[i] != 0)
			check[table[i]] ++;
	}
	for(int i = 0; i < tableSz; i++){
		if (check[i] < 1)
			printf("check[%d] = %d; too small!!!\n", i, check[i]);
	   	if (check[i] > mapNum)
		    printf("check[%d] = %d; too large!!!\n", i, check[i]);
		assert(check[i] >= 1);
		assert(check[i] <= mapNum);
	}
	printf("CHECK OK!\n");
	#endif

	if(mapNum == 1)
		for(int i = 0; i < tableSz; i++)
			S->invTable[table[i]] = i;
	else{
		if(S->factor != NULL)
			free(S->factor);
		S->factor = (double*)calloc(mapNum * tableSz, sizeof(double));
		memcpy(S->factor, factor, sizeof(double)*(mapNum * tableSz));
		if(S->invFactor != NULL)
			free(S->invFactor);
		S->invFactor = (double*)calloc(mapNum * tableSz, sizeof(double));
        for (int i=0; i<tableSz*mapNum; i++)
			S->invFactor[i] = 0;

		int count[tableSz];
		memset(count, 0, tableSz * sizeof(int));
		int elemNum = tableSz * mapNum;
		for(int i = 0; i < elemNum; i++)
			if(factor[i] != 0){
				S->invTable[table[i] * 2 + count[table[i]]] = (i / mapNum);
				S->invFactor[table[i] * 2 + count[table[i]]] = factor[i];
				count[table[i]]++;
			}
	}
}*/

//One Symmetry S may be used by many tensors,
//Be careful to do recycleSym
void recycleSym(Symmetry *S){
	free(S->table);
	free(S->invTable);
	if(S->mapNum > 1){
		free(S->factor);
		free(S->invFactor);
	}
	(S->group).~vector();
	free(S);
}
