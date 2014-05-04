#include <uni10/data-structure/Bond.h>
namespace uni10{
//namespace datatype{

Bond::Bond(bondType _type, int dim) : m_type(_type){
	Qnum q0(0);
	std::vector<Qnum> qnums(dim, q0);
	setting(qnums);
}
Bond::Bond(bondType _type, const std::vector<Qnum>& qnums) : m_type(_type){
	setting(qnums);
}

Bond::Bond(const Bond& _b):m_type(_b.m_type), m_dim(_b.m_dim), Qnums(_b.Qnums), Qdegs(_b.Qdegs), offsets(_b.offsets){
}
bondType Bond::type()const{
	return m_type;
}
int Bond::dim()const{
	return m_dim;
}
void Bond::assign(bondType _type, int dim){
	m_type = _type;
	Qnums.clear();
	Qdegs.clear();
	offsets.clear();
	Qnum q0(0);
	std::vector<Qnum> qnums(dim, q0);
	setting(qnums);
}

void Bond::assign(bondType _type, const std::vector<Qnum>& qnums){
	m_type = _type;
	Qnums.clear();
	Qdegs.clear();
	offsets.clear();
	setting(qnums);
}

void Bond::setting(const std::vector<Qnum>& qnums){
	assert(qnums.size() > 0);
	std::map<Qnum, bool> mark;
	int cnt = 0;
	m_dim = 0;
	for(int i = 0; i < qnums.size(); i++){
		if(i == 0 || !(qnums[i] == qnums[i - 1])){
			Qnums.push_back(qnums[i]);
			Qdegs.push_back(1);
			offsets.push_back(m_dim);
			cnt++;
		}
		else
			Qdegs[cnt - 1]++;
		m_dim++;
	}
}

Bond::~Bond(){}

std::ostream& operator<< (std::ostream& os, const Bond& b){
	if(b.m_type == BD_IN)
		os<<"IN : ";
	else
		os<<"OUT: ";
	for(int i = 0; i < b.Qnums.size(); i++)
		os << b.Qnums[i] << "|" << b.Qdegs[i]<<", ";
	os<<"Dim = "<< b.m_dim << std::endl;
	return os;
}

std::map<Qnum, int> Bond::degeneracy()const{
	std::map<Qnum, int>hst;
	for(int i = 0; i < Qnums.size(); i++){
		if(hst.find(Qnums[i]) == hst.end())
			hst[Qnums[i]] = Qdegs[i];
		else
			hst[Qnums[i]] += Qdegs[i];
	}
	return hst;
}
std::vector<Qnum> Bond::Qlist()const{
	std::vector<Qnum>list(m_dim);
	int cnt = 0;
	for(int q = 0; q < Qnums.size(); q++)
		for(int d = 0; d < Qdegs[q]; d++){
			list[cnt] = Qnums[q];
			cnt++;
		}
	return list;
}
bool operator== (const Bond& b1, const Bond& b2){
	return (b1.m_type == b2.m_type) && (b1.Qnums == b2.Qnums) && (b1.Qdegs == b2.Qdegs);
}
void Bond::change(bondType tp){
	if(m_type != tp){
		for(int q = 0; q < Qnums.size(); q++)
			Qnums[q] = -Qnums[q];
		m_type = tp;
	}
}
Bond& Bond::combine(Bond bd){
	bd.change(m_type);
	std::vector<Qnum> qnums;
	std::vector<int> qdegs;
	offsets.clear();
	m_dim = 0;
	Qnum qnum;
	int qdim;
	int cnt = 0;
	for(int q = 0; q < Qnums.size(); q++)
		for(int d = 0; d < Qdegs[q]; d++){
			for(int qq = 0; qq < bd.Qnums.size(); qq++){
				qnum = Qnums[q] * bd.Qnums[qq];
				qdim = bd.Qdegs[qq];
				if(qnums.size() == 0 || !(qnum == qnums[cnt - 1])){
					qnums.push_back(qnum);
					qdegs.push_back(qdim);
					offsets.push_back(m_dim);
					cnt++;
				}
				else
					qdegs[cnt - 1] += qdim;
				m_dim += qdim;
			}
		}
	Qnums = qnums;
	Qdegs = qdegs;
	return *this;
}

Bond combine(bondType tp, const std::vector<Bond>& bds){
	assert(bds.size() > 0);
	if(bds.size() > 1){
	int bd_num = bds.size();
	Bond outBond1 = bds[bd_num - 1];
	Bond outBond2 = bds[bd_num - 2];
	int b = 0;
	outBond2.change(tp);
	outBond2.combine(outBond1);
	for(b = 0; b < bd_num - 2; b++){
		if(b % 2 == 0){
			outBond1 = bds[bd_num - 3 - b];
			outBond1.change(tp);
			outBond1.combine(outBond2);
		}
		else{
			outBond2 = bds[bd_num - 3 - b];
			outBond2.change(tp);
			outBond2.combine(outBond1);
		}
	}
	if(b % 2 == 0)
		return outBond2;
	else
		return outBond1;	
	}
	else{
		return bds[0];
	}
}
Bond combine(const std::vector<Bond>& bds){
	return combine(bds[0].m_type, bds);
}
};	/* namespace uni10 */
