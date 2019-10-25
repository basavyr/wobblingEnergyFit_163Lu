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

ofstream gout("../output/Lu163.out");

typedef struct parameterSet
{
    double RMS;
    double GM;
    double I0;
    double V;
    const double V_left = 0.1;
    const double V_right = 10.0;
    const int I0_left = 40;
    const int I0_right = 60;
    const int GM_left = 0;
    const int GM_right = 60;
} params;

void showResults(params p)
{
    cout << "Results for Lu163 parameter set are:\n";
    cout << "I0   V    GAMMA \n";
    cout << p.I0 << " " << p.V << " " << p.GM << "\n";
    cout << "The RMS of the energy is = " << p.RMS << "\n";
    gout << "Results for Lu163 parameter set are:\n";
    gout << "I0   V    GAMMA \n";
    gout << p.I0 << " " << p.V << " " << p.GM << "\n";
    gout << "The RMS of the energy is = " << p.RMS << "\n";
}

EnergyCalculus Lu163;

void calculateMinimumValue(params &results)
{
    double minValue = 98765432123456789.0;
    double temp;
    for (int izero = results.I0_left; izero <= results.I0_right; ++izero)
    {
        for (int gamma = results.GM_left; gamma <= results.GM_right; ++gamma)
        {
            for (double particlePotential = results.V_left; particlePotential <= results.I0_right; particlePotential += 0.1)
            {
                temp = Lu163.rootMeanSquare((double)izero, (double)gamma * Lu163.PI / 180.0, particlePotential);
                if (temp < minValue)
                {
                    minValue = temp;
                    gout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    cout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    results.I0 = izero;
                    results.GM = gamma;
                    results.V = particlePotential;
                    results.RMS = minValue;
                }
            }
        }
    }
    Lu163.showFrequencies(results.I0, results.GM * Lu163.PI / 180.0, results.V);
}

//for reversed OMEGA ordering
void calculateMinimumValueReversed(params &results)
{
    //ofstream gout("../output/Lu163.out");
    double minValue = 98765432123456789.0;
    double temp;
    for (int izero = results.I0_left; izero <= results.I0_right; ++izero)
    {
        for (int gamma = results.GM_left; gamma <= results.GM_right; ++gamma)
        {
            for (double particlePotential = results.V_left; particlePotential <= results.I0_right; particlePotential += 0.1)
            {
                temp = Lu163.rootMeanSquareReversed((double)izero, (double)gamma * Lu163.PI / 180.0, particlePotential);
                if (temp < minValue)
                {
                    minValue = temp;
                    gout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    cout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    results.I0 = izero;
                    results.GM = gamma;
                    results.V = particlePotential;
                    results.RMS = minValue;
                }
            }
        }
    }
    Lu163.showFrequencies(results.I0, results.GM * Lu163.PI / 180.0, results.V);
}

int main()
{
    const char *appname = "Wobbling fit 2 params";
    runApp(appname);
    //sysInfo();
    params results;
    auto startTime = chrono::high_resolution_clock::now();
    /*   gout << "*** PRC ORDERING FOR OMEGA ****"
         << "\n";
    cout << "*** PRC ORDERING FOR OMEGA ****"
         << "\n";
    calculateMinimumValue(results);
    showResults(results);
    gout << "*** JPG ORDERING FOR OMEGA ****"
         << "\n";
    cout << "*** JPG ORDERING FOR OMEGA ****"
         << "\n";
    calculateMinimumValueReversed(results);
    showResults(results);
    Lu163.showEnergies(results.I0, results.GM * Lu163.PI / 180.0, results.V);
   */
    auto endTime = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    // energyDebug();
    cout << "RMS value parameter set calculation took " << (double)execTime / 1000.0 << " seconds\n";
    showDate();

    return 0;
}
