#include "twoParamsEnergies.h"
#include "expData.h"
#include "energyExpressions.h"

#include <cmath>

ExperimentalData Lu_163_Exp;
EnergyCalculus formulas_Lu163;

double TwoParams::band1EnergyShifted(double spin, double izero, double particlePotential, double shift)
{
    double eZero, energy;

    eZero = formulas_Lu163.minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (formulas_Lu163.Omega1Reversed(spinZero, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spinZero, oddSpin, izero, gamma, particlePotential));

    energy = formulas_Lu163.minHamiltonian(spin, oddSpin, izero, gamma, particlePotential) + 0.5 * (formulas_Lu163.Omega1Reversed(spin, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spin, oddSpin, izero, gamma, particlePotential));

    return energy - eZero;
}

double TwoParams::band2EnergyShifted(double spin, double izero, double particlePotential, double shift)
{
    double eZero, energy;

    eZero = formulas_Lu163.minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (formulas_Lu163.Omega1Reversed(spinZero, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spinZero, oddSpin, izero, gamma, particlePotential));

    energy = formulas_Lu163.minHamiltonian(spin - 1.0, oddSpin, izero, gamma, particlePotential) + 0.5 * (3.0 * formulas_Lu163.Omega1Reversed(spin - 1.0, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spin - 1.0, oddSpin, izero, gamma, particlePotential));

    return energy - eZero + shift;
}

double TwoParams::band3EnergyShifted(double spin, double izero, double particlePotential, double shift)
{
    double eZero, energy;

    eZero = formulas_Lu163.minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (formulas_Lu163.Omega1Reversed(spinZero, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spinZero, oddSpin, izero, gamma, particlePotential));

    energy = formulas_Lu163.minHamiltonian(spin - 2.0, oddSpin, izero, gamma, particlePotential) + 0.5 * (5.0 * formulas_Lu163.Omega1Reversed(spin - 2.0, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spin - 2.0, oddSpin, izero, gamma, particlePotential));

    return energy - eZero + shift;
}

double TwoParams::band4EnergyShifted(double spin, double izero, double particlePotential, double shift)
{
    double eZero, energy;

    eZero = formulas_Lu163.minHamiltonian(spinZero, oddSpin, izero, gamma, particlePotential) + 0.5 * (formulas_Lu163.Omega1Reversed(spinZero, oddSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spinZero, oddSpin, izero, gamma, particlePotential));

    energy = formulas_Lu163.minHamiltonian(spin - 3.0, secondSpin, izero, gamma, particlePotential) + 0.5 * (7.0 * formulas_Lu163.Omega1Reversed(spin - 3.0, secondSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spin, secondSpin, izero, gamma, particlePotential));

    return energy - eZero + singleParticleEnergy;
}

double TwoParams::shiftedSquaredSum(double izero, double particlePotential, double shift)
{
    double sum = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    bool ok = true;
    int count = 0;

    for (int i = 0; i < Lu_163_Exp.dim1 && ok; ++i)
    {
        double aux = [i, izero, particlePotential, shift](ExperimentalData exp, TwoParams th) {
            return pow(exp.spin1Exp[i] - th.band1EnergyShifted(exp.spin1Exp[i], izero, particlePotential, shift), 2.0);
        };
        if (isnan(aux))
        {
            ok = 0;
            break;
        }
        if (!isnan(aux))
        {
            sum1 += aux;
            count++;
        }
    }

    for (int i = 0; i < Lu_163_Exp.dim2 && ok; ++i)
    {
        double aux = [i, izero, particlePotential, shift](ExperimentalData exp, TwoParams th) {
            return pow(exp.spin2Exp[i] - th.band2EnergyShifted(exp.spin2Exp[i], izero, particlePotential, shift), 2.0);
        };
        if (isnan(aux))
        {
            ok = 0;
            break;
        }
        if (!isnan(aux))
        {
            sum2 += aux;
            count++;
        }
    }
    for (int i = 0; i < Lu_163_Exp.dim3 && ok; ++i)
    {
        double aux = [i, izero, particlePotential, shift](ExperimentalData exp, TwoParams th) {
            return pow(exp.spin3Exp[i] - th.band3EnergyShifted(exp.spin3Exp[i], izero, particlePotential, shift), 2.0);
        };
        if (isnan(aux))
        {
            ok = 0;
            break;
        }
        if (!isnan(aux))
        {
            sum3 += aux;
            count++;
        }
    }
    for (int i = 0; i < Lu_163_Exp.dim4 && ok; ++i)
    {
        auto aux = [i, izero, particlePotential, shift](ExperimentalData exp, TwoParams th)  {
            return (double)pow(exp.spin4Exp[i] - th.band4EnergyShifted(exp.spin4Exp[i], izero, particlePotential, shift), 2.0);
        };
        if (isnan(aux))
        {
            ok = 0;
            break;
        }
        if (!isnan(aux))
        {
            sum4 += aux;
            count++;
        }
    }
    sum = sum1 + sum2 + sum3 + sum4;
    int dim = (Lu_163_Exp.dim1 + Lu_163_Exp.dim2 + Lu_163_Exp.dim3 + Lu_163_Exp.dim4);
    if (sum && count == dim)
        return (double)sqrt(sum / dim);
    return 98765432109.0;
}
