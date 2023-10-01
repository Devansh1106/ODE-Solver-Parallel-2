#include<cmath>
#include<string>
#include<iostream>
#include <fstream>
#include <cstdio>
#include "Master.hpp"
using namespace std;

// double model_prob_1_rhs(double x){return -1.0;}
double model_prob_2_rhs(double x){return 34*sin(x);}

int main (int argc, char** argv)
{
    // SecondOrderOde ode_mp1 (1.0,0.0,0.0,model_prob_1_rhs,0.0,1.0);
    // BoundaryConditions bc_mp1;
    // bc_mp1.SetLhsDirichletBc(0.0);
    // bc_mp1.SetRhsDirichletBc(0.0);

    // BvpOde bvpode_mp1(&ode_mp1,&bc_mp1, 1001);
    // bvpode_mp1.SetFilename("model_prob_sol_1.dat");
    // bvpode_mp1.Solve();
    // Solver(argc, argv);
    
    FILE* file = fopen("filename.txt","w");
    static int flag_1 = 0;
    for (int j = 1; j < 6; j++)
    {
        BoundaryConditions bc_mp2;
        SecondOrderOde ode_mp2(1.0,3.0,-4.0,model_prob_2_rhs,0.0,M_PI);
        BvpOde bvpode_mp2(&ode_mp2,&bc_mp2,1001);
        bvpode_mp2.Rhs_gen(1001);
        fprintf(file,"model_prob_sol%d.txt\n", j);
        bc_mp2.SetLhsNeumannBc(-5.0 * j);
        bc_mp2.SetRhsDirichletBc(4.0 * j);
        if (flag_1 == 0)
        {
            bvpode_mp2.Solve();
            flag_1 = 1;         
        }

        bvpode_mp2.ApplyBoundaryConditionsRhsVector();
        bvpode_mp2.mpLinearSystem->Mumps_rhs_gen(j);
        bvpode_mp2.Clean_Rhs_gen();
        bvpode_mp2.Clean_Linearsystem();
    }
    rewind(file);
    Solver(argc, argv, file, 1001);
    return 0;
}