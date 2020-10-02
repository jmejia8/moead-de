#include <string.h>
#include <time.h>
#include "common/GlobalVariable.cpp"
// #include "ALG_EMO_MOEAD.h"
#include "moead-de.h"
void testF(vector<double> &x, vector<double> &f, const unsigned int nx)
{
	f = std::vector<double>(2, 0);

	unsigned int j;
	double sum1,g;

	sum1 = 0.0;
	for(j = 1; j <= nx-1; j++)
	{
		sum1 += x[j];
	}
	g=1+9*sum1/(nx-1);
	f[0] = x[0];
	f[1] = g*(1-sqrt(x[0]/g));
	// printf("%g, %g\n", f[0], f[1]);
}

#include "common/GlobalObject.cpp"
// #include "common/GlobalVariable.cpp"
#include "common/IndividualBase.cpp"
#include "common/SubProblemBase.cpp"
#include "common/TestInstance.cpp"
#include "common/UtilityToolBox.cpp"
#include "moead-de.cpp"




void ResetRandSeed();

int main(int argc, char const *argv[]) {

	int   NumberOfVariables;
	int   NumberOfObjectives;
	int   NumberOfFuncEvals;

	int total_run, numOfInstance;
    std::ifstream readf("Instance.txt");
	readf>>numOfInstance;
	readf>>total_run;

	char  alg_name[1024];

    //sprintf(alg_name,"MOEAD");
	sprintf(alg_name,"MOEAD-DE");

	rnd_uni_seed = 123;

	for(int inst=1; inst<=numOfInstance; inst++)
	{
		readf>>strTestInstance;
		readf>>NumberOfVariables;
		readf>>NumberOfObjectives;
		readf>>NumberOfFuncEvals;
		readf>>strCrossOverType;

		printf("\n\n -- Instance: %s, %d Variables %d  Objectives \n\n", strTestInstance, NumberOfVariables, NumberOfObjectives);

		clock_t start, temp, finish;
		double  duration, last = 0;
		start = clock();



		printf("a1\n");
		std::fstream fout;
		char logFilename[1024];
		printf("a1\n");
		sprintf(logFilename, "saving/LOG_.dat");
		printf("b1\n");
		fout.open(logFilename,std::ios::out);

		for(int run=1; run<=total_run; run++)
		{
			printf("a1\n");
			ResetRandSeed();


      if(!strcmp(alg_name,"MOEAD-DE"))
			{
			    CALG_EMO_MOEAD_DE MOEAD_DE;
					printf("start\n");
			    MOEAD_DE.Execute(run, testF,NumberOfVariables, NumberOfObjectives, NumberOfFuncEvals);
					printf("end 2\n");
			}


			temp = clock();
			duration = (double)(temp - start) / CLOCKS_PER_SEC;
			printf("time %f\n", duration - last);
			fout<<duration - last<<" ";
			last = duration;
			if(run%10==0) fout<<"\n";
			break;
		}

		fout<<"\n\n";  	finish = clock();  	fout.close();
		break;
	}

	return 0;
}

void ResetRandSeed()
{

	rnd_uni_seed = (rnd_uni_seed + 23)%1377;
	rnd_uni_init = -(long)rnd_uni_seed;

}
