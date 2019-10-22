#ifndef ENERGYEXPRESSIONS_H
#define ENEGRYEXPRESSIONS_H

//Calculates the Energy levels of each TSD band.
//Parameters: I0, V, Gamma.
class EnergyCalculus
{
private:
    const double PI = 3.141592654;
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
    double minHamiltonian(double spin, double oddSpin, double izero, double gamma, double particlePotential);

public:
    double Omega1(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    double Omega2(double spin, double oddSpin, double izero, double gamma, double particlePotential);
    double energyTSD1(double spin, double izero, double gamma, double particlePotential);
    double energyTSD2(double spin, double izero, double gamma, double particlePotential);
    double energyTSD3(double spin, double izero, double gamma, double particlePotential);
    double energyTSD4(double spin, double izero, double gamma, double particlePotential);

public:
    double squaredSum1(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum2(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum3(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double squaredSum4(double izero, double gamma, double particlePotential, const double energy[], const double spin[], int size);
    double rootMeanSquare(double izero, double gamma, double particlePotential);
};

#endif