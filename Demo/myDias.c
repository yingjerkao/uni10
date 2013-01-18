/*Function API*/
void AscendC(Tensor* W1, Tensor* W1T, Tensor* U, Tensor* H, Tensor* UT, Tensor* W2, Tensor* W2T, Tensor* Tout);
void uni_trans(Tensor* U, Tensor* UD, Tensor* H, Tensor* Tout);
/*Function Definition*/
void uni_trans(Tensor* U, Tensor* UD, Tensor* H, Tensor* Tout){
	vector<Tensor*>Tlist;
	Tlist.push_back(U);
	Tlist.push_back(UD);
	Tlist.push_back(H);
	operate(DIAG_UT, Tlist, Tout);
}
void AscendC(Tensor* W1, Tensor* W1T, Tensor* U, Tensor* H, Tensor* UT, Tensor* W2, Tensor* W2T, Tensor* Tout){
	vector<Tensor*>Tlist;
	Tlist.push_back(W1);
	Tlist.push_back(W1T);
	Tlist.push_back(U);
	Tlist.push_back(H);
	Tlist.push_back(UT);
	Tlist.push_back(W2);
	Tlist.push_back(W2T);
	operate(DIAG_AC, Tlist, Tout);
}
