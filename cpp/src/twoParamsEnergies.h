#ifndef TWOPARAMSENERGIES_H
#define TWOPARAMSENERGIES_H

class TwoParams
{
private:
    const double beta = 0.38;
    const double gamma = 17 * 3.141592654 / 180.0;
    const double spinZero = 6.5;
    const double oddSpin = 6.5;
    const double secondSpin = 6.5;
    double singleParticleEnergy = -0.334;

public:
    double band1EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band2EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band3EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band4EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double shiftedSquaredSum(double izero, double particlePotential, double shift);
    double rootMeanSquaredError(const double exp[],const double th[], int size, double P);
};

#endif // TWOPARAMSENERGIES_H
