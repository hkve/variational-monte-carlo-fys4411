#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/anharmonicoscillator.h"

#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "Solvers/metropolis.h"
#include "Solvers/metropolishastings.h"
#include "WaveFunctions/interactinggaussian.h"
#include "particle.h"
#include "sampler.h"
#include "system.h"

using namespace std;

int main(int argv, char **argc)
{
    // Seed for the random number generator
    int seed = 2023;

    // Set default paramters
    unsigned int numberOfDimensions = 3;
    unsigned int numberOfParticles = 10;
    unsigned int numberOfMetropolisSteps = (unsigned int)pow(2, 20);
    unsigned int numberOfEquilibrationSteps = (unsigned int)pow(2, 20);
    double omega = 1.0;                                 // Oscillator frequency.
    double alpha = omega / 2.0;                         // Variational parameter. If using gradient descent, this is the initial guess.
    double beta = 2.82843;                              // Beta will be 1 unless we are using the interacting wave function. Then it is 2.82843.
    double stepLength = 0.1;                            // Metropolis step length.
    double epsilon = 0.05;                              // Tolerance for gradient descent.
    double lr = 0.1;                                    // Learning rate for gradient descent.
    double interactionTerm = 0.0043 * sqrt(1. / omega); // Interection constant a. Notice that hbar = m = 1.
    bool importanceSampling = false;
    bool gradientDescent = 1;
    bool analytical = true;
    bool detailed = false;
    double D = 0.5;
    string filename = "";

    // If no arguments are given, show usage.
    if (argv == 1)
    {
        cout << "Hello! Usage:" << endl;
        cout << "./vmc #dims #particles #log10(metropolis-steps) "
                "#log10(equilibriation-steps) alpha stepLength "
                "importanceSampling? analytical? gradientDescent? detailed? filename"
             << endl;
        cout << "#dims, int: Number of dimensions" << endl;
        cout << "#particles, int: Number of particles" << endl;
        cout << "#log2(metropolis steps), int/double: log2 of number of steps, i.e. 6 gives 2^6 steps" << endl;
        cout << "#log2(@-steps), int/double: log2 of number of equilibriation steps, i.e. 6 gives 2^6 steps" << endl;
        cout << "alpha, double: WF parameter for simple gaussian. Analytical sol alpha = omega/2" << endl;
        cout << "beta, double: WF parameter for interacting gaussian. For harmonic trap, use beta=1, else beta=2.82843" << endl;
        cout << "stepLenght, double: How far should I move a particle at each MC cycle?" << endl;
        cout << "Importantce sampling?, bool: If the Metropolis Hasting algorithm is used. Then stepLength serves as Delta t" << endl;
        cout << "analytical?, bool: If the analytical expression should be used. Defaults to true" << endl;
        cout << "gradientDescent?, bool: If the gradient descent algorithm should be used. Defaults to true" << endl;
        cout << "filename, string: If the results should be dumped to a file, give the file name. If none is given, a simple print is performed." << endl;
        cout << "detailed?, bool: Spits at detail information. Defaults to " << endl;
        return 0;
    }

    // Check how many arguments are given and overwrite defaults. Works serially,
    // meaning if 4 parameters are given the first 4 paramters will be
    // overwritten, the rest will be defaults.
    if (argv >= 2)
        numberOfDimensions = (unsigned int)atoi(argc[1]);
    if (argv >= 3)
        numberOfParticles = (unsigned int)atoi(argc[2]);
    if (argv >= 4)
        numberOfMetropolisSteps = (unsigned int)pow(2, atof(argc[3]));
    if (argv >= 5)
        numberOfEquilibrationSteps = (unsigned int)pow(2, atof(argc[4]));
    if (argv >= 6)
        alpha = (double)atof(argc[5]);
    if (argv >= 7)
        beta = (double)atof(argc[6]);
    if (argv >= 8)
        stepLength = (double)atof(argc[7]);
    if (argv >= 9)
        importanceSampling = (bool)atoi(argc[8]);
    if (argv >= 10)
        analytical = (bool)atoi(argc[9]);
    if (argv >= 11)
        gradientDescent = (bool)atoi(argc[10]);
    if (argv >= 12)
        filename = argc[11];
    if (argv >= 13)
        detailed = (bool)atoi(argc[12]);

    // The random engine can also be built without a seed
    auto rng = std::make_unique<Random>(seed);

    // Initialize particles
    auto particles = setupRandomUniformInitialState(
        omega, numberOfDimensions, numberOfParticles, *rng, interactionTerm);

    // Construct a unique pointer to a new System
    std::unique_ptr<class Hamiltonian> hamiltonian;

    if (beta != 1.0)
    {
        double gamma = beta;
        hamiltonian = std::make_unique<AnharmonicOscillator>(gamma);
    }
    else
    {
        hamiltonian = std::make_unique<HarmonicOscillator>(omega);
    }

    // Initialise Interacting Gaussian by default
    std::unique_ptr<class WaveFunction> wavefunction = std::make_unique<InteractingGaussian>(
        alpha,
        beta,
        interactionTerm,
        numberOfParticles); // Empty wavefunction pointer, since it uses "alpha" in its
                            // constructor (can only be moved once).

    // Empty solver pointer, since it uses "rng" in its constructor (can only be
    // moved once).
    std::unique_ptr<class MonteCarlo> solver;

    // Set what solver to use, pass on rng and additional parameters
    if (importanceSampling)
    {
        solver = std::make_unique<MetropolisHastings>(std::move(rng), stepLength, D);
    }
    else
    {
        solver = std::make_unique<Metropolis>(std::move(rng));
    }

    // Create system pointer, passing in all classes.
    auto system = std::make_unique<System>(
        // Construct unique_ptr to Hamiltonian
        std::move(hamiltonian),
        // Construct unique_ptr to wave function
        std::move(wavefunction),
        // Construct unique_ptr to solver, and move rng
        std::move(solver),
        // Move the vector of particles to system
        std::move(particles));

    // Run steps to equilibrate particles

    if (detailed)
    {
        system->saveSamples(filename + "_blocking.dat", 0);
    }
    system->runEquilibrationSteps(stepLength, numberOfEquilibrationSteps);

    std::unique_ptr<class Sampler> sampler;
    // Run the Metropolis algorithm
    if (!gradientDescent)
    {
        sampler =
            system->runMetropolisSteps(stepLength, numberOfMetropolisSteps);
        // Output information from the simulation, either as file or print
        sampler->output(*system, filename + ".txt", omega, analytical, importanceSampling);
    }
    else
    {
        sampler = system->optimizeMetropolis(
            *system, filename, stepLength, numberOfMetropolisSteps, numberOfEquilibrationSteps, epsilon, lr);
        // Output information from the simulation, either as file or print
        sampler->output(*system, filename + ".txt", omega, analytical, importanceSampling);
    }

    if (detailed)
    {
        system->saveFinalState(filename + "_Rs.txt");
    }
    return 0;
}