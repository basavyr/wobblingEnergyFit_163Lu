#include "../include/twoParamsEnergies.h"
#include "../include/expData.h"
#include "../include/energyExpressions.h"

#include <cmath>
#include <fstream>
std::ofstream mathematicaPlots("../output/mathematicaPlots.out");

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

    energy = formulas_Lu163.minHamiltonian(spin - 3.0, secondSpin, izero, gamma, particlePotential) + 0.5 * (7.0 * formulas_Lu163.Omega1Reversed(spin - 3.0, secondSpin, izero, gamma, particlePotential) + formulas_Lu163.Omega2Reversed(spin - 3.0, secondSpin, izero, gamma, particlePotential));

    return energy - eZero + singleParticleEnergy;
}

double TwoParams::shiftedSquaredSum(double izero, double particlePotential, double shift)
{
    double sum = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    bool ok = true;
    int count = 0;
    double temp;

    //the term with index K in the RMS total sum
    auto rms = [](double exp, double th) {
        return std::pow(exp - th, 2.0);
    };

    for (int i = 0; i < Lu_163_Exp.dim1 && ok; ++i)
    {
        temp = rms(Lu_163_Exp.tsd1Exp[i], band1EnergyShifted(Lu_163_Exp.spin1Exp[i], izero, particlePotential, shift));
        if (isnan(temp))
        {
            ok = 0;
            break;
        }
        if (!isnan(temp))
        {
            count++;
            sum1 += temp;
        }
    }

    for (int i = 0; i < Lu_163_Exp.dim2 && ok; ++i)
    {
        temp = rms(Lu_163_Exp.tsd2Exp[i], band2EnergyShifted(Lu_163_Exp.spin2Exp[i], izero, particlePotential, shift));
        if (isnan(temp))
        {
            ok = 0;
            break;
        }
        if (!isnan(temp))
        {
            count++;
            sum2 += temp;
        }
    }
    for (int i = 0; i < Lu_163_Exp.dim3 && ok; ++i)
    {
        temp = rms(Lu_163_Exp.tsd3Exp[i], band3EnergyShifted(Lu_163_Exp.spin3Exp[i], izero, particlePotential, shift));
        if (isnan(temp))
        {
            ok = 0;
            break;
        }
        if (!isnan(temp))
        {
            count++;
            sum3 += temp;
        }
    }
    for (int i = 0; i < Lu_163_Exp.dim4 && ok; ++i)
    {
        temp = rms(Lu_163_Exp.tsd4Exp[i], band4EnergyShifted(Lu_163_Exp.spin4Exp[i], izero, particlePotential, shift));
        if (isnan(temp))
        {
            ok = 0;
            break;
        }
        if (!isnan(temp))
        {
            count++;
            sum4 += temp;
        }
    }
    sum = sum1 + sum2 + sum3 + sum4;
    int dim = (Lu_163_Exp.dim1 + Lu_163_Exp.dim2 + Lu_163_Exp.dim3 + Lu_163_Exp.dim4);
    if (sum && count == dim)
        return (double)sqrt(sum / dim);
    return 98765432109.0;
}

void TwoParams::plotDataMathematica(double izero, double particlePotential, double shift)
{
    mathematicaPlots << "tsd1Exp={ ";
    for (int i = 0; i < Lu_163_Exp.dim1; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin1Exp[i] << "," << Lu_163_Exp.tsd1Exp[i] << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "tsd1Th={ ";
    for (int i = 0; i < Lu_163_Exp.dim1; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin1Exp[i] << "," << band1EnergyShifted(Lu_163_Exp.spin1Exp[i], izero, particlePotential, shift) << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "\n";
    mathematicaPlots << "tsd2Exp={ ";
    for (int i = 0; i < Lu_163_Exp.dim2; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin2Exp[i] << "," << Lu_163_Exp.tsd2Exp[i] << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "tsd2Th={ ";
    for (int i = 0; i < Lu_163_Exp.dim2; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin2Exp[i] << "," << band2EnergyShifted(Lu_163_Exp.spin2Exp[i], izero, particlePotential, shift) << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "\n";
    mathematicaPlots << "tsd3Exp={ ";
    for (int i = 0; i < Lu_163_Exp.dim3; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin3Exp[i] << "," << Lu_163_Exp.tsd3Exp[i] << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "tsd3Th={ ";
    for (int i = 0; i < Lu_163_Exp.dim3; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin3Exp[i] << "," << band3EnergyShifted(Lu_163_Exp.spin3Exp[i], izero, particlePotential, shift) << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "\n";
    mathematicaPlots << "tsd4Exp={ ";
    for (int i = 0; i < Lu_163_Exp.dim4; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin4Exp[i] << "," << Lu_163_Exp.tsd4Exp[i] << " } ,";
    mathematicaPlots << "};"
                     << "\n";

    mathematicaPlots << "tsd4Th={ ";
    for (int i = 0; i < Lu_163_Exp.dim4; ++i)
        mathematicaPlots << " { " << Lu_163_Exp.spin4Exp[i] << "," << band4EnergyShifted(Lu_163_Exp.spin4Exp[i], izero, particlePotential, shift) << " } ,";
    mathematicaPlots << "};"
                     << "\n";
}