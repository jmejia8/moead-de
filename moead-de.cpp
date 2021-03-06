#include "moead-de.h"


CALG_EMO_MOEAD_DE::CALG_EMO_MOEAD_DE(void)
{

}


CALG_EMO_MOEAD_DE::~CALG_EMO_MOEAD_DE(void)
{

}


void CALG_EMO_MOEAD_DE::Execute(int run_id, void (*f)(std::vector<double>&, std::vector<double>&, unsigned int),
                                int NumberOfVariables,
                                int NumberOfObjectives,
                                int NumberOfFuncEvals)
{
    this->NumberOfVariables = NumberOfVariables;
    this->NumberOfObjectives = NumberOfObjectives;
    this->NumberOfFuncEvals = NumberOfFuncEvals;

    printf("InitializeParameter\n");

    this->InitializeParameter();
    printf("InitializePopulation\n");

    this->InitializePopulation(f);

    printf("InitializeNeighborhood\n");
    this->InitializeNeighborhood();

    int gen = 1;
    while (true) {
        gen++;

        this->EvolvePopulation(f);

        if(gen%25==0){
            printf("Instance: %s RUN: %d  GEN = %d \n", strTestInstance, run_id, gen);
        }
        if(IsTerminated()){
            break;
        }
    }


    // this->SavePopulation(run_id);

    m_PopulationSOP.clear();
    v_ReferencePoint.clear();

}



void CALG_EMO_MOEAD_DE::InitializeNeighborhood()
{
    vector<double> v_dist   = vector<double>(s_PopulationSize, 0);
    vector<int>    v_indx   = vector<int>(s_PopulationSize, 0);

    unsigned int i, j, k;
    for(i=0; i<s_PopulationSize; i++) {
        for(j=0; j<s_PopulationSize; j++) {
            v_dist[j] = UtilityToolBox.DistanceVectorNorm2(m_PopulationSOP[i].v_Weight_Vector,
                                                      				m_PopulationSOP[j].v_Weight_Vector);
            v_indx[j] = j;
        }

        UtilityToolBox.Minfastsort(v_dist, v_indx, s_PopulationSize, s_NeighborhoodSize);


        for(k=0; k<s_NeighborhoodSize; k++) {
            m_PopulationSOP[i].v_Neighbor_Index.push_back(v_indx[k]);  // save the indexes into neighborhood
        }
    }
    v_dist.clear();
    v_indx.clear();
}


void CALG_EMO_MOEAD_DE::InitializeParameter()
{
char filename[1024];

    sprintf(filename,"MOEAD-DE.txt");

    char str_temp[1024];
    std::ifstream readf(filename);

    int  flag = 0, i = 0;
    while(!readf.eof()){
        readf>>str_temp;

        if(strcmp(str_temp, strTestInstance)==0) {
            readf>>s_PopulationSize;
            readf>>s_NeighborhoodSize;
            readf>>s_LocalMatingRatio;
            readf>>s_ReplacementLimit;
            readf>>s_Scaling_Factor;
            flag = 1;
            //printf("%d %d", s_PopulationSize, s_NeighborhoodSize); getchar();
            break;
        }
    }

    readf.close();
}



void CALG_EMO_MOEAD_DE::UpdateReference(vector<double> &obj_vect)
{
    for(unsigned n=0; n<this->NumberOfObjectives; n++) {
        if(obj_vect[n]<v_ReferencePoint[n]) {
            v_ReferencePoint[n] = obj_vect[n];
        }
    }
}


void CALG_EMO_MOEAD_DE::InitializePopulation(void (*f)(std::vector<double>&, std::vector<double>&, unsigned int))
{
    unsigned i, j;

    s_Fevals_Count = 0;

    v_ReferencePoint = vector<double>(this->NumberOfObjectives, 1.0e+30);


    char filename1[1024];
    sprintf(filename1,"settings/W%dD_%d.dat", this->NumberOfObjectives, s_PopulationSize);
    std::ifstream readf(filename1);


    for(i=0; i<s_PopulationSize; i++) {
        CSubProblemBase SP(this->NumberOfVariables, this->NumberOfObjectives);

        SP.m_BestIndividual.Randomize();
        SP.m_BestIndividual.Evaluate(f);

        s_Fevals_Count++;

        UpdateReference(SP.m_BestIndividual.f_obj);    // update reference point

        for(j=0; j<this->NumberOfObjectives; j++) {
            readf>>SP.v_Weight_Vector[j];
        }

        m_PopulationSOP.push_back(SP);
    }

    readf.close();
}

void CALG_EMO_MOEAD_DE::UpdateProblem(CIndividualBase &child,
       unsigned sp_id)
{

    double f1, f2;

    int    id1 = sp_id, id2;

    vector<int> neighbor_list;
    // local update
    if(b_Selection_Local) {
        for(int i=0; i<s_NeighborhoodSize; i++) {
            neighbor_list.push_back(m_PopulationSOP[id1].v_Neighbor_Index[i]);
        }
    } else {
        // global update
        for(int i=0; i<s_PopulationSize; i++) {
            neighbor_list.push_back(i);
        }
    }


    vector<int> order(std::vector<int>(neighbor_list.size(),0));
    UtilityToolBox.RandomPermutation(order, 0);

    int count = 0;
    for(int i=0; i<neighbor_list.size(); i++) {
        id2 = neighbor_list[order[i]];
        f1 = UtilityToolBox.ScalarizingFunction(m_PopulationSOP[id2].m_BestIndividual.f_obj,
                  			                        m_PopulationSOP[id2].v_Weight_Vector,
        v_ReferencePoint,
        1);


        f2 = UtilityToolBox.ScalarizingFunction(child.f_obj,
                                        			m_PopulationSOP[id2].v_Weight_Vector,
        v_ReferencePoint,
        1);

        if(f2<f1) {
            m_PopulationSOP[id2].m_BestIndividual = child;
            count++;
        }

        if(count>=s_ReplacementLimit) {
            break;
        }
    }

    neighbor_list.clear();
    order.clear();
}




void CALG_EMO_MOEAD_DE::SelectMatingPool(vector<unsigned> &pool,
      unsigned sp_id,
      unsigned selected_size)
{
    unsigned  id, p;
    while(pool.size()<selected_size) {
        id      = (unsigned int) (s_NeighborhoodSize*UtilityToolBox.Get_Random_Number());
        p       = m_PopulationSOP[sp_id].v_Neighbor_Index[id];

        bool flag = true;
        for(unsigned i=0; i<pool.size(); i++) {
            // parent is in the list
            if(pool[i]==p) {
                flag = false;
                break;
            }
        }

        if(flag){
            pool.push_back(p);
        }
    }
}


void CALG_EMO_MOEAD_DE::EvolvePopulation(void (*f)(std::vector<double>&, std::vector<double>&, unsigned int))
{

    if(IsTerminated()) return;

    std::vector<int> order(std::vector<int>(s_PopulationSize,0));
    UtilityToolBox.RandomPermutation(order, 0);


    CIndividualBase child(this->NumberOfVariables, this->NumberOfObjectives);
    int p1, p2, p3;
    vector<unsigned> mating_pool;

    for(unsigned int s=0; s<s_PopulationSize; s++) {
        // printf("ind %d\n", s);
        unsigned int id_c = order[s];

        if(UtilityToolBox.Get_Random_Number()<s_LocalMatingRatio) {
            SelectMatingPool(mating_pool, id_c, 3);
            p1 = mating_pool[0]; p2 = mating_pool[1]; p3 = mating_pool[2]; mating_pool.clear();
            b_Selection_Local = true;
        } else {
            p1 = int(s_PopulationSize*UtilityToolBox.Get_Random_Number());
            p2 = int(s_PopulationSize*UtilityToolBox.Get_Random_Number());
            p3 = int(s_PopulationSize*UtilityToolBox.Get_Random_Number());
            b_Selection_Local = false;
        }

        /*
        UtilityToolBox.SimulatedBinaryCrossover(m_PopulationSOP[p1].m_BestIndividual.x_var,
                                            m_PopulationSOP[p2].m_BestIndividual.x_var,
        child.x_var);
        //*/

        //*
        UtilityToolBox.DifferentialEvolution(m_PopulationSOP[p1].m_BestIndividual.x_var,
                                         m_PopulationSOP[p1].m_BestIndividual.x_var,
         m_PopulationSOP[p1].m_BestIndividual.x_var,
         child.x_var,
         s_Scaling_Factor);
          //*/

        UtilityToolBox.PolynomialMutation(child.x_var, 1.0/this->NumberOfVariables);

        child.Evaluate(f);   s_Fevals_Count++;


        UpdateReference(child.f_obj);

        UpdateProblem(child, id_c);

        if(IsTerminated()) break;
    }
}


bool CALG_EMO_MOEAD_DE::IsTerminated()
{
    if(s_Fevals_Count>=this->NumberOfFuncEvals){
        return true;
    } else {
        return false;
    }

}


void CALG_EMO_MOEAD_DE::SaveObjSpace(char saveFilename[1024])
{
    std::fstream fout;
    fout.open(saveFilename,std::ios::out);
    for(unsigned n=0; n<s_PopulationSize; n++) {
        for(unsigned int k=0; k<this->NumberOfObjectives; k++) {
            fout<<m_PopulationSOP[n].m_BestIndividual.f_obj[k]<<"  ";
        }
        fout<<"\n";
    }

    fout.close();
}



void CALG_EMO_MOEAD_DE::SaveVarSpace(char saveFilename[1024])
{
    std::fstream fout;
    fout.open(saveFilename,std::ios::out);
    for(unsigned n=0; n<s_PopulationSize; n++) {
        for(unsigned k=0; k<this->NumberOfVariables; k++) {
            fout<<m_PopulationSOP[n].m_BestIndividual.x_var[k]<<"  ";
        }
        fout<<"\n";
    }
    fout.close();
}

void CALG_EMO_MOEAD_DE::SavePopulation(int run_id)
{
    char filename[1024];
    sprintf(filename,"saving/POF_%s_RUN%d.dat",strTestInstance, run_id);
    SaveObjSpace(filename);

    sprintf(filename,"saving/POS_%s_RUN%d.dat",strTestInstance, run_id);
    SaveVarSpace(filename);
}
