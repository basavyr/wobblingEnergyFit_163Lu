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
    showDate();
    return 0;
}
