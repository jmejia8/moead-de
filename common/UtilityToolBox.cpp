// UtilityTool.cpp: implementation of the CUtilityToolBox class.
//
//////////////////////////////////////////////////////////////////////

#include "UtilityToolBox.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CUtilityToolBox::CUtilityToolBox()
{

}

CUtilityToolBox::~CUtilityToolBox()
{

}


double CUtilityToolBox::DistanceVectorPBI(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
	double sum = 0.;
	double lambda = InnerProduct(vec1, vec2)/InnerProduct(vec2, vec2);
	for(int i=0; i<dim; i++)
		sum += (vec1[i] - lambda*vec2[i])*(vec1[i] - lambda*vec2[i]);
	return sum;
}

void CUtilityToolBox::NormalizeVector(vector<double> &vect)
{
	double sum = 0;
	for(int i=0; i<vect.size(); i++)
	{
	   sum += vect[i];
	}

	for(int j=0; j<vect.size(); j++)
	{
	   vect[j] = vect[j]/sum;
	}
}

double CUtilityToolBox::InnerProduct(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=vec1[n]*vec2[n];
	return sum;
}

double CUtilityToolBox::DistanceVectorNorm2(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double CUtilityToolBox::VectorNorm2(vector <double> &vec1)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=vec1[n]*vec1[n];
	return sqrt(sum);
}

double CUtilityToolBox::DistanceVectorNorm1(vector<double> &vec1, 
											vector<double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=fabs(vec1[n] - vec2[n]);
	return sum;
}


double CUtilityToolBox::ScalarizingFunction(vector<double> &y_obj, 
											vector<double> &namda, 
											vector<double> &referencepoint, 
											int s_type)
{

	// Chebycheff Scalarizing Function

	unsigned n, nobj = y_obj.size(); 

	double max_fun = -1.0e+30, diff, feval, scalar_obj;     

	if(s_type==1)
	{
		for(n=0; n<nobj; n++)
		{
			diff = fabs(y_obj[n] - referencepoint[n] + 0.01);
			if(namda[n]==0) 
				feval = 0.001*diff;
			else
			    feval = namda[n]*diff;

			if(feval>max_fun){
				max_fun = feval;
			}
		}
		scalar_obj = max_fun;
	}


	if(s_type==2)
	{
		double alpha = 0.1;

		for(n=0; n<nobj; n++)
		{
			diff  = y_obj[n] - referencepoint[n];
			feval = diff*(1.0/nobj + alpha*namda[n]);
			if(feval>max_fun){
				max_fun = feval;
			}
		}
		scalar_obj = max_fun;
	}

	if(s_type==3)   // PBI
	{
		vector<double> vect_normal = vector<double>(nobj,0);
		double namda_norm = this->VectorNorm2(namda);
		for(int n=0; n<nobj; n++)
		{
			vect_normal[n] = namda[n]/namda_norm;
		}			

		double d1_proj = abs(this->InnerProduct(y_obj, vect_normal) - this->InnerProduct(referencepoint, vect_normal));
		double temp = d1_proj*d1_proj-2*d1_proj*(this->InnerProduct(y_obj,vect_normal) - this->InnerProduct(referencepoint,vect_normal))
			          + (this->InnerProduct(y_obj, y_obj) + this->InnerProduct(referencepoint, referencepoint) - 2*this->InnerProduct(y_obj, referencepoint));				
		double d2_pred  = sqrt(temp);
		scalar_obj = d1_proj + 20*d2_pred;
	}
	return scalar_obj;
}


void CUtilityToolBox::RandomPermutation(vector<int> &permutation, int type)
{
	if(type==0)
	{
	    for(unsigned k=0; k<permutation.size(); k++) permutation[k] = k;
	}
	random_shuffle(permutation.begin(), permutation.end());

}


vector<int> CUtilityToolBox::IndexOfMinimumInt(vector<int> &vec)
{
	vector<int> I;
	I.clear();
	int id = -1;
	int min_value = 1e9; // 整型最大为2^32
	for(int i=0; i<vec.size(); i++)
	{
		if(vec[i]<min_value)
		{
			min_value = vec[i];
			id = i;
		}
	}
	for(int i=0; i<vec.size(); i++)
		if(vec[i]==min_value)
			I.push_back(i);
	return I;
}

int CUtilityToolBox::IndexOfMinimumDouble(vector<double> &vec)
{
	int id = 0;
	double minvalue = 1.0e30;
	for(int i=0; i<vec.size(); i++)
	{
		if(vec[i]<minvalue)
		{
			minvalue = vec[i];
			id = i;
		}
	}
	return id;
}

void CUtilityToolBox::Minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}
void CUtilityToolBox::Maxfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]<x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}



void CUtilityToolBox::PolynomialMutation(vector<double>& x_var, double rate)
{
    double rand1, rand2, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy, lowBound = 0, uppBound = 1;
	double eta_m = 20;

	int nvar = x_var.size();

    for(unsigned int j=0; j<nvar; j++)
    {
		rand1 = Get_Random_Number();
        if (rand1<=rate)
        {
            y  = x_var[j];
            yl = lowBound;
            yu = uppBound;
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rand2 = Get_Random_Number();
            mut_pow = 1.0/(eta_m+1.0);
            if (rand2 <= 0.5)
            {
                xy  = 1.0 - delta1;
                val = 2.0*rand2+(1.0 - 2.0*rand2)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy  = 1.0 - delta2;
                val = 2.0*(1.0-rand2)+2.0*(rand2 - 0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            
			x_var[j] = y;
        }
    }
    return;
}

void CUtilityToolBox::SimulatedBinaryCrossover(vector<double> &x_var1, vector<double> &x_var2, vector<double> &x_var3)
{
    double rnd;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
	double eta_c = 20;
	double lowBound = 0, 
		   uppBound = 1;
	int    nvar = x_var1.size();
    if (Get_Random_Number()<= 1.0) 
    {
        for (int i=0; i<nvar; i++)
        {
            if (Get_Random_Number()<=0.5 )
            {
                if (fabs(x_var1[i] - x_var2[i]) > EPS)
                {
                    if (x_var1[i] < x_var2[i])
                    {
                        y1 = x_var1[i];
                        y2 = x_var2[i];
                    }
                    else
                    {
                        y1 = x_var2[i];
                        y2 = x_var1[i];
                    }
                    yl  = lowBound;
                    yu  = uppBound;
                    rnd = Get_Random_Number();
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rnd <= (1.0/alpha))
                    {
                        betaq = pow ((rnd*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rnd*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rnd <= (1.0/alpha))
                    {
                        betaq = pow ((rnd*alpha),(1.0/(eta_c+1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0/(2.0 - rnd*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    //if(rand()<=0.5)
					if(Get_Random_Number()<=0.5)
                    {
                        x_var3[i] = c2;
                        //child2.x_var[i] = c1;
                    }
                    else
                    {
                        x_var3[i] = c1;
                        //child2.x_var[i] = c2;
                    }
                }
                else
                {
                    x_var3[i] = x_var1[i];
                    //child2.x_var[i] = parent2.x_var[i];
                }
            }
            else
            {
                x_var3[i] = x_var1[i];
                //child2.x_var[i] = parent2.x_var[i];
            }
        }
    }
    else
    {
        for (int i=0; i<nvar; i++)
        {
            x_var3[i] = x_var1[i];
            //child2.x_var[i] = parent2.x_var[i];
        }
    }
    return;
}

void CUtilityToolBox::DifferentialEvolution(vector<double> &x_var0, vector<double> &x_var1, vector<double> &x_var2, vector<double> &x_var, double rate)
{
	double rand0 = Get_Random_Number();
	int idx_rnd  = int(rand0*x_var.size());

	double lowBound = 0, uppBound = 1;
	
	for(unsigned int n=0; n<x_var.size(); n++)
	{
	  /*Selected Two Parents*/

	  // strategy one 
	  // x_var[n] = x_var0[n] + rate*(x_var2[n] - x_var1[n]);
	  
	  //*
	  // strategy two
	  double rand1 = Get_Random_Number();

	  double CR    = 1.0;
	  if(rand1<CR||n==idx_rnd)
		  x_var[n] = x_var0[n] + rate*(x_var2[n] - x_var1[n]);
	  else
		  x_var[n] = x_var0[n];
	  //*/


	  // handle the boundary voilation
	  if(x_var[n]<lowBound)
	  {
		  double rand2 = Get_Random_Number();
		  x_var[n]     = lowBound + rand2*(x_var0[n] - lowBound);
	  }
	  if(x_var[n]>uppBound)
	  { 
		  double rand3 = Get_Random_Number();
		  x_var[n]     = uppBound - rand3*(uppBound - x_var0[n]);
	  }
	}
}

int CUtilityToolBox::GetWeightNumber(int nobj, int H)
{
    int a = H + nobj - 1, b = nobj - 1, prod1 = 1, prod2 = 1;

	
	for(int i=1; i<=b; i++)      prod1*=(a-i+1);
	for(int j=1; j<=b; j++)      prod2*=j;


	return int(1.0*prod1/prod2);
}


double CUtilityToolBox::MaxElementInVector(vector<double> &vect)
{
	double max_e = vect[0];
	for(unsigned n=1; n<vect.size(); n++)
	{
	    if(vect[n]>max_e)
		{
		    max_e = vect[n];
		}
	}
	return max_e;
}


double CUtilityToolBox::Get_Random_Number()
{
    return Rnd_Uni(&rnd_uni_init);
}


double CUtilityToolBox::Rnd_Uni(long *idum)
{

    long j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0)
	{
	if (-(*idum) < 1) *idum=1;
	else *idum = -(*idum);
	idum2=(*idum);
	for (j=NTAB+7;j>=0;j--)
	{
	  k=(*idum)/IQ1;
	  *idum=IA1*(*idum-k*IQ1)-k*IR1;
	  if (*idum < 0) *idum += IM1;
	  if (j < NTAB) iv[j] = *idum;
	}
	iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;

}