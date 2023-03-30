#include "system.h"

#include <cassert>
#include <iostream>
#include <memory>
#include <iomanip>
#include <string>

#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Solvers/montecarlo.h"
#include "WaveFunctions/wavefunction.h"
#include "particle.h"
#include "sampler.h"
#include <fstream>

using namespace std;
System::System(std::unique_ptr<class Hamiltonian> hamiltonian,
               std::unique_ptr<class WaveFunction> waveFunction,
               std::unique_ptr<class MonteCarlo> solver,
               std::vector<std::unique_ptr<class Particle>> particles)
{
  m_numberOfParticles = particles.size();
  m_numberOfDimensions = particles[0]->getNumberOfDimensions();
  m_hamiltonian = std::move(hamiltonian);
  m_waveFunction = std::move(waveFunction);
  m_solver = std::move(solver);
  m_particles = std::move(particles);
}

unsigned int System::runEquilibrationSteps(
    double stepLength, unsigned int numberOfEquilibrationSteps)
{
  unsigned int acceptedSteps = 0;

  for (unsigned int i = 0; i < numberOfEquilibrationSteps; i++)
  {
    acceptedSteps += m_solver->step(stepLength, *m_waveFunction, m_particles);
  }

  return acceptedSteps;
}

std::unique_ptr<class Sampler> System::runMetropolisSteps(
    double stepLength, unsigned int numberOfMetropolisSteps)
{
  auto sampler =
      std::make_unique<Sampler>(m_numberOfParticles, m_numberOfDimensions,
                                stepLength, numberOfMetropolisSteps);

  if (m_saveSamples)
    sampler->openSaveSample(m_saveSamplesFilename);

  for (unsigned int i = 0; i < numberOfMetropolisSteps; i++)
  {
    /* Call solver method to do a single Monte-Carlo step.
     */
    bool acceptedStep =
        m_solver->step(stepLength, *m_waveFunction, m_particles);

    // compute local energy
    sampler->sample(acceptedStep, this);

    if (m_saveSamples)
      sampler->saveSample(i);
  }

  sampler->computeAverages();
  if (m_saveSamples)
    sampler->closeSaveSample();

  return sampler;
}

std::unique_ptr<class Sampler> System::optimizeMetropolis(
    System &system, std::string filename, double stepLength, unsigned int numberOfMetropolisSteps, unsigned int numberOfEquilibrationSteps,
    double epsilon, double learningRate)
{
  double gradient = 1;
  int epoch = 0;
  double alpha_0 = getWaveFunctionParameters()[0];
  double alpha, beta;
  double beta_0 = getWaveFunctionParameters()[1];
  auto sampler =
      std::make_unique<Sampler>(m_numberOfParticles, m_numberOfDimensions,
                                stepLength, numberOfMetropolisSteps);

  // run equilibration steps and store positions into vector
  runEquilibrationSteps(stepLength, numberOfEquilibrationSteps);

  for (unsigned int i = 0; i < m_numberOfParticles; i++)
  {
    m_particles[i]->saveEquilibrationPosition(); // by doind this, we just need to do equilibriation once in the GD
  }

  while (std::abs(gradient) > epsilon)
  {

    gradient = 0;
    /*Positions are reset to what they were after equilibration, but the parameters of the wave function should be what they were at the END of last epoch*/

    // reset position and quantum force (quantum force is reset automatically to 0 at the beggining of the metroplis step)

    for (unsigned int i = 0; i < m_numberOfParticles; i++)
    {
      m_particles[i]->resetEquilibrationPosition();
    }

    sampler = system.runMetropolisSteps(stepLength, numberOfMetropolisSteps);

    alpha = getWaveFunctionParameters()[0];
    beta = getWaveFunctionParameters()[1]; // this might be useless since we are not asked to optimize beta
    std::vector<double> parameters = getWaveFunctionParameters();
    sampler->writeGradientSearchToFile(system, filename, alpha_0, epoch, alpha, beta_0, beta);

    std::vector<double> m_energyDerivative = sampler->getEnergyDerivative();

    int n_params = m_waveFunction->getNumberOfParameters();

    for (int i = 0; i < n_params - 1; i++) // we do not want to optimize beta because minima is given to us already. But if we want, the functionality is there!
    {
      parameters[i] -= learningRate * m_energyDerivative[i];
      std::cout << "parameters post update: " << parameters[i] << "\n";
      std::cout << "m_energyDerivative: " << m_energyDerivative[i] << "\n";
      gradient += m_energyDerivative[i];
    }
    epoch++;
    // set new wave function parameters
    m_waveFunction->setParameters(parameters);
  }

  alpha = getWaveFunctionParameters()[0];
  beta = getWaveFunctionParameters()[1];
  //  get the last epoch values
  sampler->writeGradientSearchToFile(system, filename, alpha_0, epoch, alpha, beta_0, beta);

  return sampler;
}

double System::computeLocalEnergy()
{
  // Helper function
  return m_hamiltonian->computeLocalEnergy(*m_waveFunction, m_particles);
}

const std::vector<double> &System::getWaveFunctionParameters()
{
  // Helper function
  return m_waveFunction->getParameters();
}

void System::setWaveFunction(std::unique_ptr<class WaveFunction> waveFunction)
{
  m_waveFunction = std::move(waveFunction);
}

void System::setSolver(std::unique_ptr<class MonteCarlo> solver)
{
  m_solver = std::move(solver);
}

double System::computeParamDerivative(int paramIndex)
{
  // Helper function
  return m_waveFunction->computeParamDerivative(m_particles, paramIndex);
}

void System::saveSamples(std::string filename, int skip)
{
  // Tells system to save local energy estimates during run
  m_saveSamples = true;
  m_saveSamplesFilename = filename;

  // Due to powers of two being just a single bit 1, use bitwise AND to check if skip is a power of 2.
  bool isPow2 = false;
  if (skip == 0)
    isPow2 = true;
  else
    isPow2 = skip > 0 && !(skip & (skip - 1));

  assert(isPow2);
  m_skip = skip;
}

int System::getSkip()
{
  return m_skip;
}

void System::saveFinalState(std::string filename)
{
    std::ofstream file(filename, std::ios::out | std::ios::trunc);

    int w = 20;
    file << setw(w) << "x" << setw(w) << "y" << setw(w) << "z\n";

    for(int i = 0; i < m_numberOfParticles; i++)
    {
        auto r = m_particles.at(i)->getPosition();
        file << setw(w) << r.at(0) << setw(w) << r.at(1) << setw(w) << r.at(2) << "\n";
    } 
    file.close();
}