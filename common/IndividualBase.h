// Individual.h: interface for the CIndividualBase class.
//
//////////////////////////////////////////////////////////////////////
#pragma once

#include "GlobalVariable.h"
#include "GlobalObject.h"
#include <vector>

using namespace std;

class CIndividualBase
{
public:
	CIndividualBase();
	CIndividualBase(int, int);
	virtual ~CIndividualBase();

	vector <double> x_var;
	vector <double> f_obj;
	vector <double> f_normal;


	unsigned int    rank;
	unsigned int    count;
	unsigned int    type;
	double          density;
	int NumberOfVariables;
	int NumberOfObjectives;

	void initialize(int, int);
	void   Randomize();
	void   Evaluate(void (*f)(std::vector<double>&, std::vector<double>&, unsigned int));
	void   Show(int type);

    bool   operator<(const CIndividualBase &ind2);
	bool   operator<<(const CIndividualBase &ind2);
    bool   operator==(const CIndividualBase &ind2);
    void   operator=(const CIndividualBase &ind2);

};
