#include <iostream>
#include <uni10.hpp>

int main(){
	uni10::Qnum q1(1);	
	uni10::Qnum q0(0);	
	uni10::Qnum q_1(-1);	
	// Create an array of Qnums for the states of a bond.
	std::vector<uni10::Qnum> qnums;
	qnums.push_back(q1);
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);

	// Constrcut Bond with Qnum array
	uni10::Bond bd(uni10::BD_IN, qnums);
	// Print out a Bond
	std::cout<<"Bond bd: \n"<<bd<<std::endl;
	std::cout<<"Bond type: "<<bd.type()<<"(IN)"<<", Bond dimension: "<<bd.dim()<<std::endl<<std::endl;

	// List the degeneracy of states
	std::cout<<"Degeneracies: "<<std::endl;
	std::map<uni10::Qnum, int> degs = bd.degeneracy();
	for(std::map<uni10::Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it)
		std::cout<<it->first<<": "<<it->second<<std::endl;
	std::cout<<std::endl;

	std::vector<uni10::Qnum> qlist = bd.Qlist();
	std::cout<<"Qnum list: "<<std::endl;
	for(int i = 0; i < qlist.size(); i++)
		std::cout<<qlist[i]<<", ";
	std::cout<<std::endl<<std::endl;

	// Change bond type
	bd.change(uni10::BD_OUT);
	std::cout<<"bd changes to BD_OUT:\n"<<bd<<std::endl;

	bd.change(uni10::BD_IN);

	// Combine bond
	qnums.clear();
	qnums.push_back(q1);
	qnums.push_back(q0);
	qnums.push_back(q0);
	qnums.push_back(q_1);
	uni10::Bond bd2(uni10::BD_IN, qnums);
	std::cout<<"Bond bd2: \n"<<bd2<<std::endl;

	// bd.combine(bd2);
	std::cout<<"bd2.combine(bd): \n"<<bd2.combine(bd)<<std::endl;

	std::cout<<"Degeneracies of bd2 after combining bd: "<<std::endl;
	degs = bd2.degeneracy();
	for(std::map<uni10::Qnum,int>::const_iterator it=degs.begin(); it!=degs.end(); ++it)
		std::cout<<it->first<<": "<<it->second<<std::endl;
	std::cout<<std::endl;

	return 0;
}


