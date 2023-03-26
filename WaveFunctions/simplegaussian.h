#pragma once

#include <memory>

#include "wavefunction.h"

class SimpleGaussian : public WaveFunction
{
public:
    SimpleGaussian(double alpha, double beta);
    double evaluate(std::vector<std::unique_ptr<class Particle>> &particles);
    double evaluate_w(int proposed_particle_idx, class Particle &proposed_particle, class Particle &old_particle, std::vector<std::unique_ptr<class Particle>> &particles);
    double computeParamDerivative(std::vector<std::unique_ptr<class Particle>> &particles, int parameterIndex);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles);
    void quantumForce(std::vector<std::unique_ptr<class Particle>> &particles, Particle &particle, std::vector<double> &force);
    void setParameters(std::vector<double> parameters);
};

class SimpleGaussianNumerical : public SimpleGaussian
{
public:
    SimpleGaussianNumerical(double alpha, double beta, double dx);
    double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles);
    double evaluateSingleParticle(class Particle &particle);

private:
    double m_dx;
};