#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>

#include "common/UtilityToolBox.h"
#include "common/GlobalVariable.h"
#include "common/IndividualBase.h"
#include "common/SubProblemBase.h"



class CALG_EMO_MOEAD_DE
{
public:
	CALG_EMO_MOEAD_DE(void);
	~CALG_EMO_MOEAD_DE(void);

	void Execute(int run_id);

	void InitializeNeighborhood();
	void InitializePopulation();
	void InitializeParameter();



	void UpdateReference(vector<double> &obj_vect);
    void UpdateProblem(CIndividualBase &child, unsigned sp_id);

    void SelectMatingPool(vector<unsigned> &pool, unsigned sp_id, unsigned selected_size);
	void EvolvePopulation();
	bool IsTerminated();


	void SaveObjSpace(char saveFilename[1024]);
	void SaveVarSpace(char saveFilename[1024]);
	void SavePopulation(int run_id);


public:

    vector <CSubProblemBase> m_PopulationSOP;
	vector <double>      v_ReferencePoint;

	unsigned int     s_PopulationSize;
    unsigned int	 s_NeighborhoodSize;
	unsigned int     s_ReplacementLimit;
	double           s_LocalMatingRatio;
	double           s_Scaling_Factor;

	int              s_Fevals_Count;
	bool             b_Selection_Local;

};
