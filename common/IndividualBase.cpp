// Individual.cpp: implementation of the CIndividualBase class.
//
//////////////////////////////////////////////////////////////////////

#include "IndividualBase.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CIndividualBase::CIndividualBase()
{
    x_var    = vector<double>(NumberOfVariables, 0);
    f_obj    = vector<double>(NumberOfObjectives, 0);
    f_normal = vector<double>(NumberOfObjectives, 0);
}

CIndividualBase::~CIndividualBase()
{

}

void CIndividualBase::Randomize()
{
    int lowBound = 0, uppBound = 1;

    unsigned int n;

    for(n=0; n<NumberOfVariables; n++)
    {
        x_var[n] = lowBound + UtilityToolBox.Get_Random_Number()*(uppBound - lowBound);
    }
}

void CIndividualBase::Evaluate()
{


    if(!strcmp("ZDT1", strTestInstance))   TestInstance.ZDT1(x_var, f_obj, x_var.size());
    if(!strcmp("ZDT2", strTestInstance))   TestInstance.ZDT2(x_var, f_obj, x_var.size());
    if(!strcmp("ZDT3", strTestInstance))   TestInstance.ZDT3(x_var, f_obj, x_var.size());

    if(!strcmp("DTLZ1", strTestInstance))  TestInstance.DTLZ1(x_var, f_obj, x_var.size());
    if(!strcmp("DTLZ2", strTestInstance))  TestInstance.DTLZ2(x_var, f_obj, x_var.size());


    if(!strcmp("LZ1", strTestInstance))  TestInstance.TEC09_LZ1(x_var, f_obj, x_var.size());
    if(!strcmp("LZ2", strTestInstance))  TestInstance.TEC09_LZ2(x_var, f_obj, x_var.size());
    if(!strcmp("LZ3", strTestInstance))  TestInstance.TEC09_LZ3(x_var, f_obj, x_var.size());
    if(!strcmp("LZ4", strTestInstance))  TestInstance.TEC09_LZ4(x_var, f_obj, x_var.size());
    if(!strcmp("LZ5", strTestInstance))  TestInstance.TEC09_LZ5(x_var, f_obj, x_var.size());
    if(!strcmp("LZ6", strTestInstance))  TestInstance.TEC09_LZ6(x_var, f_obj, x_var.size());
    if(!strcmp("LZ7", strTestInstance))  TestInstance.TEC09_LZ7(x_var, f_obj, x_var.size());
    if(!strcmp("LZ8", strTestInstance))  TestInstance.TEC09_LZ8(x_var, f_obj, x_var.size());
    if(!strcmp("LZ9", strTestInstance))  TestInstance.TEC09_LZ9(x_var, f_obj, x_var.size());

}


void CIndividualBase::Show(int type)
{

    unsigned int n;
    if(type==0)
    {
        for(n=0; n<NumberOfObjectives; n++)
        printf("%f ",f_obj[n]);
        printf("\n");
    }
    else
    {
        for(n=0; n<NumberOfVariables; n++)
        printf("%f ",x_var[n]);
        printf("\n");
    }
}


void CIndividualBase::operator=(const CIndividualBase &ind2)
{
    x_var = ind2.x_var;
    f_obj = ind2.f_obj;
    f_normal = ind2.f_normal;
    rank  = ind2.rank;
    type  = ind2.type;
    count = ind2.count;
    density = ind2.density;
}

bool CIndividualBase::operator<(const CIndividualBase &ind2)
{
    bool dominated = true;
    unsigned int n;
    for(n=0; n<NumberOfObjectives; n++)
    {
        if(ind2.f_obj[n]<f_obj[n]) return false;
    }
    if(ind2.f_obj==f_obj) return false;
    return dominated;
}


bool CIndividualBase::operator<<(const CIndividualBase &ind2)
{
    bool dominated = true;
    unsigned int n;
    for(n=0; n<NumberOfObjectives; n++)
    {
        if(ind2.f_obj[n]<f_obj[n]  - 0.0001) return false;
    }
    if(ind2.f_obj==f_obj) return false;
    return dominated;
}

bool CIndividualBase::operator==(const CIndividualBase &ind2)
{
    if(ind2.f_obj==f_obj) return true;
    else return false;
}
