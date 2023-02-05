#include <memory>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "system.h"
#include "sampler.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


Sampler::Sampler(
        unsigned int numberOfParticles,
        unsigned int numberOfDimensions,
        double stepLength,
        unsigned int numberOfMetropolisSteps)
{
    m_stepNumber = 0;
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_numberOfParticles = numberOfParticles;
    m_numberOfDimensions = numberOfDimensions;
    
    m_energy = 0;
    m_energy_variance = 0;
    m_energy_std = 0;
    
    m_cumulativeEnergy = 0;
    m_cumulativeEnergy2 = 0;
    m_numberOfAcceptedSteps = 0;
}


void Sampler::sample(bool acceptedStep, System* system) {
    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */
    auto localEnergy = system->computeLocalEnergy();
    m_cumulativeEnergy  += localEnergy;
    m_cumulativeEnergy2 += (localEnergy*localEnergy);
    m_stepNumber++;
    m_numberOfAcceptedSteps += acceptedStep;
}

void Sampler::printOutputToTerminal(System& system) {
    auto pa = system.getWaveFunctionParameters();
    auto p = pa.size();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << m_numberOfParticles << endl;
    cout << " Number of dimensions : " << m_numberOfDimensions << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(m_numberOfMetropolisSteps) << endl;
    cout << " Ratio of accepted steps: " << ((double) m_numberOfAcceptedSteps) / ((double) m_numberOfMetropolisSteps) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (unsigned int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " Energy variance : " << m_energy_variance << endl;
    cout << " Energy std : " << m_energy_std << endl;
    cout << endl;
}

void Sampler::writeOutToFile(System& system, std::string filename) {
    
}

void Sampler::computeAverages() {
    /* Compute the averages of the sampled quantities.
     */
    m_energy = m_cumulativeEnergy / m_numberOfMetropolisSteps;
     
    m_cumulativeEnergy2 /= m_numberOfMetropolisSteps;
    m_energy_variance = (m_cumulativeEnergy2 - m_energy*m_energy);
    m_energy_std = sqrt(m_energy_variance);
}
