#ifndef ENERGYEXPRESSIONS_H
#define ENEGRYEXPRESSIONS_H

//Calculates the Energy levels of each TSD band.
//Parameters: I0, V, Gamma.
class EnergyCalculus
{
private:
    const double BETA = 0.38;
    double a1, a2, a3;
    double inertiaMoment(double izero, double gamma, int k);
    double inertiaFactor1(double izero, double gamma);
    double inertiaFactor2(double izero, double gamma);
    double inertiaFactor3(double izero, double gamma);
    double subTerm1(double spin, double oddSpin, double izero, double gamma);
    double subTerm2(double spin, double oddSpin, double izero, double gamma);
    double subTerm3(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    double subTerm4(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    double bigTerm1(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    double bigTerm2(double spin, double oddSpin, double izero, double gamma, double particlePotential);

public:
    const double PI = 3.141592654;
    double minHamiltonian(double spin, double oddSpin, double izero, double gamma, double particlePotential);

    //order of the frequencies from PRC2017
    double Omega1(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    //order of the frequencies from PRC2017
    double Omega2(double spin, double oddSpin, double izero, double gamma, double particlePotential);

    //order of the frequencies from JPG2018
    double Omega1Reversed(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    //order of the frequencies from JPG2018
    double Omega2Reversed(double spin, double oddSpin, double izero, double gamma, double particlePotential);

    //order of the frequencies from PRC2017
    double energyTSD1(double spin, double izero, double gamma, double particlePotential);
    double energyTSD2(double spin, double izero, double gamma, double particlePotential);
    double energyTSD3(double spin, double izero, double gamma, double particlePotential);
    double energyTSD4(double spin, double izero, double gamma, double particlePotential);

    //order of the frequencies from JPG2018
    double energyTSD1Reversed(double spin, double izero, double gamma, double particlePotential);
    double energyTSD2Reversed(double spin, double izero, double gamma, double particlePotential);
    double energyTSD3Reversed(double spin, double izero, double gamma, double particlePotential);
    double energyTSD4Reversed(double spin, double izero, double gamma, double particlePotential);

    //order of the frequencies from PRC2017
public:
    double squaredSum1(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum2(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum3(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum4(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double rootMeanSquare(double izero, double gamma, double particlePotential);

    //order of the frequencies from JPG2018
public:
    double squaredSum1Reversed(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum2Reversed(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum3Reversed(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum4Reversed(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double rootMeanSquareReversed(double izero, double gamma, double particlePotential);

    //Debugging function
public:
    void squaredSumDebug1(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result);
    void squaredSumDebug2(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result);
    void squaredSumDebug3(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result);
    void squaredSumDebug4(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size, double &result);

public:
    void showFrequencies(double izero, double gamma, double particlePotential);
    void showEnergies(double izero, double gamma, double particlePotential);
};

#endif