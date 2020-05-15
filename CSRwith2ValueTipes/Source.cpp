#include <iostream>
#include"sstream"
#include<memory>
#include"cg.h"
#include"CSRClasses.h"

using namespace std;

int main()
{
	cout << "Enter number, indicating file's name" << endl;	
	// შეგვყავს 3, ან 3600, ან ფაილის ნომერი
	// 3, or 3600, or any number, which coincides with the name of the data file
	int k;
	cin >> k;
	stringstream converter;
	converter << k;
	string name{ converter.str() };
	name += ".txt";

	//run compromised version
	unique_ptr<simpleAbstractSM> m = make_unique<GenericCSR>(name);
	cg::multiSolve(cg::cgSymmPosDefPrecond<simpleAbstractSM>, move(m), 1e-5);
}