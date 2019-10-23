#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>
#include <fstream>

//optional OMP
#include <omp.h>

//headers
#include "expData.h"
#include "sysInfo.h"
#include "energyExpressions.h"

using namespace std;

void calculateMinimumValue()
{
    ofstream gout("../output/Lu163.out");
    double minValue = 98765432123456789.0;
    EnergyCalculus Lu163;
    double temp;
    for (int izero = 1; izero <= 200; ++izero)
    {
        for (int gamma = 0; gamma <= 60; ++gamma)
        {
            for (double particlePotential = 0.01; particlePotential <= 20; particlePotential += 0.1)
            {
                temp = Lu163.rootMeanSquare((double)izero, (double)gamma * Lu163.PI / 180.0, particlePotential);
                if (temp < minValue)
                {
                    minValue = temp;
                    gout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    cout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                }
            }
        }
    }
}

//debugging the energy formula
void energyDebug()
{
    EnergyCalculus formulas;
    ExperimentalData Lu163;
    double res;
    cout << "TSD1"
         << "\n";
    formulas.squaredSumDebug1(11, 0.12, 4, Lu163.tsd1Exp, Lu163.spin1Exp, Lu163.dim1, res);
    cout << "RMS result= " << res << "\n";
    cout << "TSD2"
         << "\n";
    formulas.squaredSumDebug2(11, 0.12, 4, Lu163.tsd2Exp, Lu163.spin2Exp, Lu163.dim2, res);
    cout << "RMS result= " << res << "\n";
    cout << "TSD3"
         << "\n";
    formulas.squaredSumDebug3(11, 60, 4, Lu163.tsd3Exp, Lu163.spin3Exp, Lu163.dim3, res);
    cout << "RMS result= " << res << "\n";
    cout << "TSD4"
         << "\n";
    formulas.squaredSumDebug4(11, 0.12, 4, Lu163.tsd4Exp, Lu163.spin4Exp, Lu163.dim4, res);
    cout << "RMS result= " << res << "\n";
}

int main()
{
    const char *appname = "Wobbling fit 2 params";
    runApp(appname);
    //sysInfo();
    calculateMinimumValue();
    // energyDebug();
    showDate();
    return 0;
}
