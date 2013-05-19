#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <string>
using namespace std;
#include "TensorLib.h"

int main(){
	string str;
	ifstream infile;
	infile.open ("AscendL");
	int lnum = 7;
	int l;
	int pos;
	int endpos;
	string tar("1234567890-");
	vector<string> names;
	vector< vector<int> > label_arr;
	for(l = 0; l < lnum; l++){
		getline(infile, str); // Saves the line in STRING.
		if(infile.eof())
			break;
		pos = str.find(":");
		names.push_back(str.substr(0, pos));
		vector<int> labels;
		while((pos = str.find_first_of(tar, pos + 1)) != string::npos){
			endpos = str.find_first_not_of(tar, pos + 1);
			string label;
			if(endpos == string::npos)
				label = str.substr(pos);
			else
				label = str.substr(pos, endpos - pos);
			char* pEnd;
			labels.push_back(strtol(label.c_str(), &pEnd, 10));
			pos = endpos;
			if(pos == string::npos)
				break;
		}
		label_arr.push_back(labels);
	}

	assert(l == lnum);	
	for(int i = 0; i < names.size(); i++){
		cout<<names[i]<<": ";
		for(int j = 0; j < label_arr[i].size(); j++)
			cout<<label_arr[i][j]<<" ";
		cout<<endl;
	}
	infile.close();
}

