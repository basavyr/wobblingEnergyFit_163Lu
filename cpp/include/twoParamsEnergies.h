#ifndef TWOPARAMSENERGIES_H
#define TWOPARAMSENERGIES_H

class TwoParams
{
private:
    const double beta = 0.38;
    const double spinZero = 6.5;
    const double oddSpin = 6.5;
    const double secondSpin = 4.5;
    double singleParticleEnergy = -0.334;

public:
    const double gamma = 17.0 * 3.141592654 / 180.0;
    double band1EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band2EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band3EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double band4EnergyShifted(double spin, double izero, double particlePotential, double shift);
    double shiftedSquaredSum(double izero, double particlePotential, double shift);

    
    //MATHEMATICA DATA FORMAT FOR PLOTTING
public:
    void plotDataMathematica(double izero, double particlePotential,double shift);
};

#endif // TWOPARAMSENERGIES_H
