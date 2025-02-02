#pragma once

#include <memory>
#include <vector>

#include "../particle.h"
#include "Math/random.h"


std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& randomEngine
            );

std::vector<std::unique_ptr<Particle>> setupRandomUniformInitialState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& randomEngine,
            double a
            );

std::vector<std::unique_ptr<Particle>> setupRandomGaussianState(
            double omega,
            unsigned int numberOfDimensions,
            unsigned int numberOfParticles,
            Random& randomEngine
            );