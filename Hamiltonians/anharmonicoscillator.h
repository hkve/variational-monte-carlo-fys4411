#pragma once
#include <memory>
#include <vector>

#include "hamiltonian.h"

class AnharmonicOscillator : public Hamiltonian
{
public:
    AnharmonicOscillator(double gamma);
    double computeLocalEnergy(
        class WaveFunction &waveFunction,
        std::vector<std::unique_ptr<class Particle>> &particles);

private:
    double m_gamma;
};
