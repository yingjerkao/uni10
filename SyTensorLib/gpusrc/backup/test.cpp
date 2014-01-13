#include <stdio.h>
#include <map>
#include <iostream>
using namespace std;

int main(){
	map<int, int> test;
	test[3] = 6;
	test[5] = 10;
	cout<<test[3]<<endl;
	cout<<test[5]<<endl;
	cout<<test[4]<<endl;

}
