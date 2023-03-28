#include <memory>
#include <cassert>
#include <iostream>

#include "anharmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

AnharmonicOscillator::AnharmonicOscillator(double gamma)
{
    m_gamma = gamma;
    assert(m_gamma == 2.82843);
}

double AnharmonicOscillator::computeLocalEnergy(
    class WaveFunction &waveFunction,
    std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* Here, you need to compute the kinetic and potential energies.
     * Access to the wave function methods can be done using the dot notation
     * for references, e.g., wavefunction.computeDoubleDerivative(particles),
     * to get the Laplacian of the wave function.
     * */
    int num_particles = particles.size();
    int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    // double psi_T = waveFunction.evaluate(particles);
    double r2_xy = 0;
    double r2_z = 0;
    double r_q = 0;

    for (int k = 0; k < num_particles; k++)
    {
        Particle &particle = *particles.at(k);
        for (int q = 0; q < numberOfDimensions - 1; q++)
        {
            r_q = particle.m_position.at(q);
            r2_xy += r_q * r_q;
        }
        r_q = particle.m_position.at(numberOfDimensions - 1);
        r2_z += r_q * r_q;
    }

    double potentialEnergy = 0.5 * (r2_xy + m_gamma * m_gamma * r2_z);

    // kinectic energy does not change with anharmonic oscillator
    double kineticEnergy = -0.5 * waveFunction.computeDoubleDerivative(particles);

    return kineticEnergy + potentialEnergy;
}
