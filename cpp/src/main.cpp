#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>

//optional OMP
#include <omp.h>

//headers
#include "expData.h"
#include "sysInfo.h"
#include "energyExpressions.h"

using namespace std;

int main()
{
    const char *appname = "Wobbling fit 2 params";
    runApp(appname);
    //sysInfo();
    EnergyCalculus Lu163;
    for (int izero = 50; izero <= 60; ++izero)
    {
        for (int gamma = 10; gamma <= 11; ++gamma)
        {
            cout << (double)izero << " " << double(gamma) << " " << Lu163.rootMeanSquare((double)izero, (double)gamma * Lu163.PI / 180.0, 0.5) << "\n";
        }
    }
    showDate();
    return 0;
}
