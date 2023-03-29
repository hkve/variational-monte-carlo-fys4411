#include <memory>
#include <iostream>
#include <cassert>
#include <cmath>

#include "initialstate.h"
#include "../particle.h"
#include "Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
        )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();
    double characteristicLength = 1/std::sqrt(omega);

    for (unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (unsigned int j=0; j < numberOfDimensions; j++) {
            double q = (rng.nextDouble()-0.5)*characteristicLength;
            position.push_back(q);
        }

        particles.push_back(std::make_unique<Particle>(position));
    }

    return particles;
}

std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng,
            double a
        )
{
    bool valid = true;
    std::vector<std::unique_ptr<Particle>> particles;
    unsigned int i, j;
    do 
    {
        valid = true;
        particles = setupRandomUniformInitialState(omega, numberOfDimensions, numberOfParticles, rng);


        for(i = 0; i < numberOfParticles; i++)
        {
            for(j = i+1; j < numberOfParticles; j++)
            {
                if( particle_r2(*particles[i], *particles[j]) < a*a )
                {

                    valid = false;
                    std::cout << "Particle distance configured to < a = " << a << ". Reconfiguring\n"; 
                    break;
                }
            }
            if(!valid)
                break;
        }
    } while(!valid);

    return particles;
}

std::vector<std::unique_ptr<Particle>> setupRandomGaussianState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& rng
            )
{
    assert(numberOfDimensions > 0 && numberOfParticles > 0);

    auto particles = std::vector<std::unique_ptr<Particle>>();
    double characteristicLength = 1/std::sqrt(omega);
   
   for(unsigned int i=0; i < numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for(unsigned int j=0; j < numberOfDimensions; j++) {
            double q = rng.nextGaussian(0.0, 1.0)*characteristicLength;
            position.push_back(q);
        }

        particles.push_back(std::make_unique<Particle>(position));
   }

   return particles;
}