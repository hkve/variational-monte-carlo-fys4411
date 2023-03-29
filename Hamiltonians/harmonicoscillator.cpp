#include <memory>
#include <cassert>
#include <iostream>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(double omega)
{
    assert(omega > 0);
    m_omega = omega;
}

double HarmonicOscillator::computeLocalEnergy(
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
    double r2_sum = 0;
    double r_q = 0;
    
    for (int k = 0; k < num_particles; k++)
    {
        Particle &particle = *particles.at(k);
        r2_sum += particle_r2(particle);
    }

    double potentialEnergy = 0.5 * r2_sum; // 0.5 * omega^2 * r^2 but omega = 1 always!!
    double kineticEnergy = -0.5 * waveFunction.computeDoubleDerivative(particles);

    return kineticEnergy + potentialEnergy;
}
