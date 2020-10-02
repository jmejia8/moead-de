#include "SubProblemBase.h"


CSubProblemBase::CSubProblemBase(int NumberOfVariables, int NumberOfObjectives)
{
		// CIndividualBase tmp(int NumberOfVariables, int NumberOfObjectives);
		m_BestIndividual.initialize(NumberOfVariables, NumberOfObjectives);// = tmp;
    v_Weight_Vector   = std::vector<double>(NumberOfObjectives, 0);
	v_Weight_Vector2   = std::vector<double>(NumberOfObjectives, 0);
}

CSubProblemBase::~CSubProblemBase()
{

}

void CSubProblemBase::Show_Weight_Vector()
{
   unsigned int n;
   for(n=0; n<v_Weight_Vector.size(); n++)
   {
       printf("%f ", v_Weight_Vector[n]);
   }
   printf("\n");
}

void CSubProblemBase::Normalize_Weight_Vector()
{
   double sum = 0;
   unsigned int s;
   for(s=0; s<v_Weight_Vector.size(); s++)
   {
       sum += v_Weight_Vector[s];
   }

   for(s=0; s<v_Weight_Vector.size(); s++)
   {
	   v_Weight_Vector[s] /= sum;
   }
}

void CSubProblemBase::operator=(const CSubProblemBase &SP)
{
    m_BestIndividual       = SP.m_BestIndividual;
	v_Neighbor_Index       = SP.v_Neighbor_Index;
	v_Weight_Vector        = SP.v_Weight_Vector;
	s_Count               = SP.s_Count;
}
