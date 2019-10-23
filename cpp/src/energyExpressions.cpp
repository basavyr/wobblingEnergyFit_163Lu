#include "energyExpressions.h"
#include "expData.h"

#include <cmath>
#include <iostream>

double EnergyCalculus::inertiaMoment(double izero, double gamma, int k)
{
    return (double)izero / (1.0 + sqrt(5.0 / (16.0 * PI)) * BETA) * (1.0 - sqrt(5.0 / (4.0 * PI)) * BETA * cos(gamma + 2.0 / 3.0 * PI * k));
}

double EnergyCalculus::inertiaFactor1(double izero, double gamma)
{
    return (double)1.0 / (2.0 * inertiaMoment(izero, gamma, 1));
}

double EnergyCalculus::inertiaFactor2(double izero, double gamma)
{
    return (double)1.0 / (2.0 * inertiaMoment(izero, gamma, 2));
}

double EnergyCalculus::inertiaFactor3(double izero, double gamma)
{
    return (double)1.0 / (2.0 * inertiaMoment(izero, gamma, 3));
}

double EnergyCalculus::subTerm1(double spin, double oddSpin, double izero, double gamma)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    return (2.0 * spin - 1.0) * (a3 - a1) + 2.0 * oddSpin * a1;
}

double EnergyCalculus::subTerm2(double spin, double oddSpin, double izero, double gamma)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    return (2.0 * spin - 1.0) * (a2 - a1) + 2.0 * oddSpin * a1;
}

double EnergyCalculus::subTerm3(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    double sing, cosg, rad3;
    rad3 = sqrt(3.0);
    sing = sin(gamma);
    cosg = cos(gamma);
    return (2.0 * oddSpin - 1.0) * (a3 - a1) + 2.0 * spin * a1 + particlePotential * (2.0 * oddSpin - 1.0) / (oddSpin * (oddSpin + 1.0)) * rad3 * (rad3 * cosg + sing);
}

double EnergyCalculus::subTerm4(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    double sing, cosg, rad3;
    rad3 = sqrt(3.0);
    sing = sin(gamma);
    cosg = cos(gamma);
    return (2.0 * oddSpin - 1.0) * (a2 - a1) + 2.0 * spin * a1 + particlePotential * (2.0 * oddSpin - 1.0) / (oddSpin * (oddSpin + 1.0)) * 2.0 * rad3 * sing;
}

double EnergyCalculus::bigTerm1(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    return -1.0 * (subTerm1(spin, oddSpin, izero, gamma) * subTerm2(spin, oddSpin, izero, gamma) + 8.0 * a2 * a3 * spin * oddSpin + subTerm3(spin, oddSpin, izero, gamma, particlePotential) * subTerm4(spin, oddSpin, izero, gamma, particlePotential));
}

double EnergyCalculus::bigTerm2(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    return (subTerm1(spin, oddSpin, izero, gamma) * subTerm3(spin, oddSpin, izero, gamma, particlePotential) - 4.0 * spin * oddSpin * a3 * a3) * (subTerm2(spin, oddSpin, izero, gamma) * subTerm4(spin, oddSpin, izero, gamma, particlePotential) - 4.0 * spin * oddSpin * a2 * a2);
}

double EnergyCalculus::Omega1(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    return sqrt(0.5 * (-bigTerm1(spin, oddSpin, izero, gamma, particlePotential) - pow(pow(bigTerm1(spin, oddSpin, izero, gamma, particlePotential), 2) - 4.0 * bigTerm2(spin, oddSpin, izero, gamma, particlePotential), 0.5)));
}
double EnergyCalculus::Omega2(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    return sqrt(0.5 * (-bigTerm1(spin, oddSpin, izero, gamma, particlePotential) + pow(pow(bigTerm1(spin, oddSpin, izero, gamma, particlePotential), 2) - 4.0 * bigTerm2(spin, oddSpin, izero, gamma, particlePotential), 0.5)));
}

double EnergyCalculus::minHamiltonian(double spin, double oddSpin, double izero, double gamma, double particlePotential)
{
    double sing;
    sing = sin(gamma + PI / 6.0);
    a1 = inertiaFactor1(izero, gamma);
    a2 = inertiaFactor2(izero, gamma);
    a3 = inertiaFactor3(izero, gamma);
    return (a2 + a3) * (spin + oddSpin) / 2.0 + a1 * pow(spin - oddSpin, 2) - particlePotential * (2.0 * oddSpin - 1.0) / (oddSpin + 1.0) * sing;
}

double EnergyCalculus::energyTSD1(double spin, double izero, double gamma, double particlePotential)
{
    double rez, eZero, spinZero, oddSpin;
    spinZero = 6.5;
    oddSpin = 6.5;

    //excitation energy of the first TSD1 band I=13/2
    eZero = minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (Omega1(spinZero, oddSpin, izero, gamma, particlePotential) + Omega2(spinZero, oddSpin, izero, gamma, particlePotential));

    //ACTUAL ENERGY OF THE BAND LEVEL I
    rez = minHamiltonian(spin, oddSpin, izero, gamma, particlePotential) + 0.5 * (Omega1(spin, oddSpin, izero, gamma, particlePotential) + Omega2(spin, oddSpin, izero, gamma, particlePotential));

    return rez - eZero;
}

double EnergyCalculus::energyTSD2(double spin, double izero, double gamma, double particlePotential)
{
    double rez, eZero, spinZero, oddSpin;
    spinZero = 6.5;
    oddSpin = 6.5;

    //excitation energy of the first TSD1 band I=13/2
    eZero = minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (Omega1(spinZero, oddSpin, izero, gamma, particlePotential) + Omega2(spinZero, oddSpin, izero, gamma, particlePotential));

    //ACTUAL ENERGY OF THE BAND LEVEL I
    rez = minHamiltonian(spin - 1.0, oddSpin, izero, gamma, particlePotential) + 0.5 * (3.0 * Omega1(spin - 1.0, oddSpin, izero, gamma, particlePotential) + Omega2(spin - 1.0, oddSpin, izero, gamma, particlePotential));

    return rez - eZero;
}

double EnergyCalculus::energyTSD3(double spin, double izero, double gamma, double particlePotential)
{
    double rez, eZero, spinZero, oddSpin;
    spinZero = 6.5;
    oddSpin = 6.5;

    //excitation energy of the first TSD1 band I=13/2
    eZero = minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (Omega1(spinZero, oddSpin, izero, gamma, particlePotential) + Omega2(spinZero, oddSpin, izero, gamma, particlePotential));

    //ACTUAL ENERGY OF THE BAND LEVEL I
    rez = minHamiltonian(spin - 2.0, oddSpin, izero, gamma, particlePotential) + 0.5 * (5.0 * Omega1(spin - 2.0, oddSpin, izero, gamma, particlePotential) + Omega2(spin - 2.0, oddSpin, izero, gamma, particlePotential));

    return rez - eZero;
}

double EnergyCalculus::energyTSD4(double spin, double izero, double gamma, double particlePotential)
{
    double rez, eZero, spinZero, oddSpin, secondSpin;
    double singleParticleEnergy = -0.334;
    spinZero = 6.5;
    oddSpin = 6.5;
    secondSpin = 4.5;

    //excitation energy of the first TSD1 band I=13/2
    eZero = minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (Omega1(spinZero, oddSpin, izero, gamma, particlePotential) + Omega2(spinZero, oddSpin, izero, gamma, particlePotential));

    //ACTUAL ENERGY OF THE BAND LEVEL I
    rez = minHamiltonian(spin - 3.0, secondSpin, izero, gamma, particlePotential) + 0.5 * (7.0 * Omega1(spin - 3.0, secondSpin, izero, gamma, particlePotential) + Omega2(spin - 3.0, secondSpin, izero, gamma, particlePotential));

    return rez - eZero + singleParticleEnergy;
}

double EnergyCalculus::squaredSum1(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD1(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            if (isnan(temp))
            {
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    if (sum && count == size)
        return (double)sum / size;
    return 9876543210;
}

double EnergyCalculus::squaredSum2(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD2(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            if (isnan(temp))
            {
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    if (sum && count == size)
        return (double)sum / size;
    return 9876543210;
}

double EnergyCalculus::squaredSum3(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD3(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            if (isnan(temp))
            {
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    if (sum && count == size)
        return (double)sum / size;
    return 9876543210;
}

double EnergyCalculus::squaredSum4(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD4(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            if (isnan(temp))
            {
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    if (sum && count == size)
        return (double)sum / size;
    return 9876543210;
}

double EnergyCalculus::rootMeanSquare(double izero, double gamma, double particlePotential)
{
    double sum1, sum2, sum3, sum4;
    ExperimentalData Lu163;
    sum1 = squaredSum1(izero, gamma, particlePotential, Lu163.tsd1Exp, Lu163.spin1Exp, Lu163.dim1);
    sum2 = squaredSum2(izero, gamma, particlePotential, Lu163.tsd2Exp, Lu163.spin2Exp, Lu163.dim2);
    sum3 = squaredSum3(izero, gamma, particlePotential, Lu163.tsd3Exp, Lu163.spin3Exp, Lu163.dim3);
    sum4 = squaredSum4(izero, gamma, particlePotential, Lu163.tsd4Exp, Lu163.spin4Exp, Lu163.dim4);
    return (double)sum1 + sum2 + sum3 + sum4;
}

//Debugging

void EnergyCalculus::squaredSumDebug1(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD1(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            std::cout << count << " " << i << " " << temp << " " << energy[i] << " " << energyTSD1(spin[i], izero, gamma, particlePotential) << "\n";
            if (isnan(temp))
            {
                std::cout << "is nan TRUE for temp = " << temp << " and index= " << i << "\n";
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                std::cout << "is nan FALSE for temp = " << temp << " and index= " << i << "\n";
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    result = 9876543210.0;
    if (sum && count == size)
    {
        std::cout << "SUM CHECKS!"
                  << "\n";
        std::cout << (double)sum / size << "\n";
        result = (double)sum / size;
    }
}

void EnergyCalculus::squaredSumDebug2(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD2(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            std::cout << count << " " << i << " " << temp << " " << energy[i] << " " << energyTSD2(spin[i], izero, gamma, particlePotential) << "\n";
            if (isnan(temp))
            {
                std::cout << "is nan TRUE for temp = " << temp << " and index= " << i << "\n";
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                std::cout << "is nan FALSE for temp = " << temp << " and index= " << i << "\n";
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    result = 9876543210.0;
    if (sum && count == size)
    {
        std::cout << "SUM CHECKS!"
                  << "\n";
        std::cout << (double)sum / size << "\n";
        result = (double)sum / size;
    }
}
void EnergyCalculus::squaredSumDebug3(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD3(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            std::cout << count << " " << i << " " << temp << " " << energy[i] << " " << energyTSD3(spin[i], izero, gamma, particlePotential) << "\n";
            if (isnan(temp))
            {
                std::cout << "is nan TRUE for temp = " << temp << " and index= " << i << "\n";
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                std::cout << "is nan FALSE for temp = " << temp << " and index= " << i << "\n";
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    result = 9876543210.0;
    if (sum && count == size)
    {
        std::cout << "SUM CHECKS!"
                  << "\n";
        std::cout << (double)sum / size << "\n";
        result = (double)sum / size;
    }
}
void EnergyCalculus::squaredSumDebug4(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result)
{
    double sum = 0.0;
    double temp;
    bool ok = true;
    int count = 0;
    while (ok)
    {
        for (int i = 0; i < size && ok; ++i)
        {
            temp = pow(energy[i] - energyTSD4(spin[i], izero, gamma, particlePotential), 2.0) / energy[i];
            std::cout << count << " " << i << " " << temp << " " << energy[i] << " " << energyTSD4(spin[i], izero, gamma, particlePotential) << "\n";
            if (isnan(temp))
            {
                std::cout << "is nan TRUE for temp = " << temp << " and index= " << i << "\n";
                ok = false;
                break;
            }
            if (!isnan(temp))
            {
                sum += temp;
                std::cout << "is nan FALSE for temp = " << temp << " and index= " << i << "\n";
                count++;
            }
        }
        if (ok)
            ok = false;
    }
    result = 9876543210.0;
    if (sum && count == size)
    {
        std::cout << "SUM CHECKS!"
                  << "\n";
        std::cout << (double)sum / size << "\n";
        result = (double)sum / size;
    }
}