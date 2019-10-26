#include "../include/expData.h"
#include <iostream>
#include <cstring>

void ExperimentalData::showData(const double arr[], const int size)
{
    std::string message = "1";
    switch (size)
    {
    case dim1:
        message = "TSD1";
        break;
    case dim2:
        message = "TSD2";
        break;
    case dim3:
        message = "TSD3";
        break;
    case dim4:
        message = "TSD4";
        break;
    default:
        break;
    }
    std::cout << "******** " << message << " *******"
              << "\n";
    for (int i = 0; i < size; ++i)
        std::cout << arr[i] << " ";
    std::cout << "\n";
}

void ExperimentalData::showBand(const double energy[], const double spin[], const int size)
{
    std::string message = "1";
    switch (size)
    {
    case dim1:
        message = "TSD1";
        break;
    case dim2:
        message = "TSD2";
        break;
    case dim3:
        message = "TSD3";
        break;
    case dim4:
        message = "TSD4";
        break;
    default:
        break;
    }
    std::cout << "******** " << message << " *******"
              << "\n";
    std::cout << "I [hbar]     "
              << "E [MeV]"
              << "\n";
    for (int i = 0; i < size; ++i)
        std::cout << spin[i] << "  " << energy[i] << "\n";
}
