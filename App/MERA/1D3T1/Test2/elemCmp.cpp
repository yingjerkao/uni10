bool elemCmp(SyTensor_t& ten1, SyTensor_t& ten2){
	int n1 = ten1.getElemNum();
	int n2 = ten2.getElemNum();
	double diff;
	if(n1 == n2){
		for(int i = 0; i < n1; i++){
			diff = fabs(ten1.elem[i] - ten2.elem[i]);
			if(diff > 1E-6)
				return false;
		}
	}
	else
		return false;
	return true;
}
