#include "expData.h"
#include <iostream>

void ExperimentalData::showData(double arr[])
{
    for (int i = 0; i <= dim1; ++i)
        std::cout << tsd1Exp[i] << " ";
    std::cout << "\n";
}

int main()
{
    ExperimentalData myClass;
    for (auto &n : myClass.spin1Exp)
        std::cout << n << " ";
    myClass.showData();
    return 0;
}
