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

int main(int argc, char const *argv[]) {

	int   NumberOfVariables = 10;
	int   NumberOfObjectives = 2;
	int   NumberOfFuncEvals = 1000;

	CALG_EMO_MOEAD_DE MOEAD_DE;
	printf("start\n");
	MOEAD_DE.Execute(1, testF,NumberOfVariables, NumberOfObjectives, NumberOfFuncEvals);

	return 0;
}
