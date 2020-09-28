// UtilityTool.h: interface for the CUtilityToolBox class.
//
//////////////////////////////////////////////////////////////////////

#pragma once


#include <vector>
#include <algorithm>
#include <math.h>
#include "GlobalVariable.h"

using namespace std;

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


class CUtilityToolBox  
{
public:
	CUtilityToolBox();
	virtual ~CUtilityToolBox();


public:
	double ScalarizingFunction(vector<double> &y_obj, 
		                      vector<double> &namda, 
							  vector<double> &referencepoint, 
							  int s_type);



    void   RandomPermutation(vector<int> &permutation, int type);

	vector<int> IndexOfMinimumInt(vector<int> &vec);
	int         IndexOfMinimumDouble(vector<double> &vec);

	void   Minfastsort(vector<double> &x, vector<int> &idx, int n, int m);
	void   Maxfastsort(vector<double> &x, vector<int> &idx, int n, int m);
	void   Minfastsort(vector<int> &x, vector<int> &idx, int n, int m);
    
	void   PolynomialMutation(vector<double>& x_var, double m_rate);
    void   SimulatedBinaryCrossover(vector<double> &x_var1, 
		                            vector<double> &x_var2, 
									vector<double> &child);
    void   DifferentialEvolution(vector<double> &x_var0, 
		                         vector<double> &x_var1, 
								 vector<double> &x_var2, 
								 vector<double> &child, 
								 double rate);
	
	
	int    GetWeightNumber(int nobj, int H);


	double DistanceVectorNorm1(vector <double> &vec1, vector <double> &vec2);	
    double DistanceVectorNorm2(vector <double> &vec1, vector <double> &vec2);
	double VectorNorm2(vector <double> &vec1);
	void   NormalizeVector(vector<double> &vect);

	double DistanceVectorPBI(vector <double> &vec1, vector <double> &vec2);
	double InnerProduct(vector <double> &vec1, vector <double> &vec2);
	
	double MaxElementInVector(vector<double> &vect);	


    double Rnd_Uni(long *idum);
	double Get_Random_Number();
};

