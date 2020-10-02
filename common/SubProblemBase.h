#pragma once

#include "GlobalVariable.h"
#include "GlobalObject.h"
#include "IndividualBase.h"
#include <vector>

using namespace std;


class CSubProblemBase
{

public:

	CIndividualBase   m_BestIndividual;
	vector<int>       v_Neighbor_Index;
	vector<double>    v_Weight_Vector;
	vector<double>    v_Weight_Vector2;
	int               s_Count;

public:

     CSubProblemBase();
     CSubProblemBase(int, int);
    ~CSubProblemBase();

	void Show_Weight_Vector();
	void Normalize_Weight_Vector();
	void operator=(const CSubProblemBase &SP);

};


class CSubProblemMAP:public CSubProblemBase
{
public:

     CSubProblemMAP();
    ~CSubProblemMAP();

public:
	vector<CIndividualBase> m_Population;

	void operator=(const CSubProblemMAP &SP);
};
