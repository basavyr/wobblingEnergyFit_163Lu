#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

//optional OMP
#include <omp.h>

//headers
#include "../include/expData.h"
#include "../include/sysInfo.h"
#include "../include/energyExpressions.h"
#include "../include/twoParamsEnergies.h"

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
    const double V_left = 0.1;
    const double V_right = 10.0;
    const int I0_left = 40;
    const int I0_right = 60;
    const int GM_left = 0;
    const int GM_right = 60;
    vector<double> shiftedRMS;
    double GM;
    double I0;
    double RMS;
    double V;
    vector<double> shiftedGM;
    vector<double> shiftedI0;
    vector<double> shiftedV;
} params;

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
                }
            }
        }
    }
    //Lu163.showFrequencies(results.I0, results.GM * Lu163.PI / 180.0, results.V);
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
                temp = Lu163.rootMeanSquareReversed((double)izero, (double)17.0 * Lu163.PI / 180.0, particlePotential);
                if (temp < minValue)
                {
                    minValue = temp;
                    gout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                    cout << izero << " " << gamma << " " << particlePotential << " " << minValue << "\n";
                }
            }
        }
    }
    //Lu163.showFrequencies(results.I0, results.GM * Lu163.PI / 180.0, results.V);
}

void shifteTwoParamsMinimum(params &results, double shift)
{
    double minValue = 98765432123456789.0;
    double temp = 0;
    TwoParams twoparams;
    for (double izero = 30; izero <= 100; izero += 1)
    {
        for (double V = 0.1; V <= 20.0; V += 0.1)
        {
            temp = twoparams.shiftedSquaredSum(izero, V, shift);
            if (temp < minValue)
            {
                minValue = temp;
                //cout << izero << " " << V << " " << minValue << "\n";
                //gout << izero << " " << V << " " << minValue << "\n";
                results.I0 = izero;
                results.V = V;
                results.RMS = minValue;
            }
        }
    }
}

ofstream shiftedOutout("../output/shiftedEnergies.out");

void showOutput(params P)
{
    cout << "I0 = " << P.I0 << "\n";
    cout << "V = " << P.V << "\n";
    cout << "RMS = " << P.RMS << "\n";
}

void generatePlotData(double izero, double particlePotential, double shift)
{
    TwoParams params;
    params.plotDataMathematica(izero, particlePotential, shift);
}

void keepRunningIt(int k)
{
    params results;
    while (k >= 0)
    {
        auto shift = 0.1 * k;
        shifteTwoParamsMinimum(results, shift);
        //showOutput(results);
        results.shiftedRMS.push_back(results.RMS);
        results.shiftedI0.push_back(results.I0);
        results.shiftedV.push_back(results.V);
        cout << "\n";
        k--;
    }
    shiftedOutout << "{ ";
    for (auto &n : results.shiftedRMS)
    {
        shiftedOutout << n << " , ";
    }
    shiftedOutout << " } \n ";
    auto minElementIndex = std::min_element(results.shiftedRMS.begin(), results.shiftedRMS.end()) - results.shiftedRMS.begin();
    cout << "Minimal set of shifted parameters is: \n";
    cout << "I0 = " << results.shiftedI0.at(minElementIndex) << "\n";
    cout << "V = " << results.shiftedV.at(minElementIndex) << "\n";
    cout << "RMS = " << results.shiftedRMS.at(minElementIndex) << "\n";
    shiftedOutout << "Minimal set of shifted parameters is: \n";
    shiftedOutout << "I0 = " << results.shiftedI0.at(minElementIndex) << "\n";
    shiftedOutout << "V = " << results.shiftedV.at(minElementIndex) << "\n";
    shiftedOutout << "RMS = " << results.shiftedRMS.at(minElementIndex) << "\n";
    generatePlotData(results.shiftedI0.at(minElementIndex), results.shiftedV.at(minElementIndex), results.shiftedRMS.at(minElementIndex));
}

int main()
{
    const char *appname = "Wobbling fit 2 params";
    cout << "\n";
    runApp(appname);
    cout << "\n";
    //sysInfo();
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
    keepRunningIt(15);
    auto endTime = chrono::high_resolution_clock::now();
    auto execTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    // energyDebug();
    cout << "\n";
    cout << "RMS value parameter set calculation took " << (double)execTime / 1000.0 << " seconds\n";
    cout << "\n";
    showDate();

    return 0;
}
