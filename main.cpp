#include "common/GlobalVariable.cpp"
// #include "ALG_EMO_MOEAD.h"
#include "moead-de.h"

#include "common/GlobalObject.cpp"
// #include "common/GlobalVariable.cpp"
#include "common/IndividualBase.cpp"
#include "common/SubProblemBase.cpp"
#include "common/TestInstance.cpp"
#include "common/UtilityToolBox.cpp"
#include "moead-de.cpp"

#include <string.h>
#include <time.h>

void ResetRandSeed();

int main(int argc, char const *argv[]) {

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

		printf("\n\n -- Instance: %s, %d Variables %d  Objectives \n\n", strTestInstance, NumberOfVariables,
						                         NumberOfObjectives);

		clock_t start, temp, finish;
		double  duration, last = 0;
		start = clock();



		std::fstream fout;
		char logFilename[1024];
		// sprintf(logFilename, "../../SAVING/%s/LOG/LOG_%s.dat", alg_name, strTestInstance);
		fout.open(logFilename,std::ios::out);

		for(int run=1; run<=total_run; run++)
		{
			ResetRandSeed();

			if(!strcmp(alg_name,"MOEAD"))
			{
				printf("nada por aqui\n");
				return 0;
			    // CALG_EMO_MOEAD MOEAD;
			    // MOEAD.Execute(run);
					//
			}

            if(!strcmp(alg_name,"MOEAD-DE"))
			{
			    CALG_EMO_MOEAD_DE MOEAD_DE;
			    MOEAD_DE.Execute(run);
			}


			temp = clock();
			duration = (double)(temp - start) / CLOCKS_PER_SEC;
			fout<<duration - last<<" ";
			last = duration;
			if(run%10==0) fout<<"\n";
		}

		fout<<"\n\n";  	finish = clock();  	fout.close();
	}

	return 0;
}

void ResetRandSeed()
{

	rnd_uni_seed = (rnd_uni_seed + 23)%1377;
	rnd_uni_init = -(long)rnd_uni_seed;

}