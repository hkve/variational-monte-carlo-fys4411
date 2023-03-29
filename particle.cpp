#include "particle.h"
#include <cassert>

Particle::Particle(const std::vector<double> &position)
{
    m_numberOfDimensions = position.size();
    m_position = position;
    m_initialPosition = position; // Save initial position. Notice that this is a copy, not a reference.
    m_EquilibrationPosition = position;
}

void Particle::adjustPosition(double change, unsigned int dimension)
{
    m_position[dimension] += change;
}
void Particle::setPosition(double new_position, unsigned int dimension)
{
    m_position[dimension] = new_position;
}
void Particle::resetPosition()
{
    m_position = m_initialPosition;
}

double particle_r2(Particle &p)
{
    /*
    Calculate r^2 for particle p
    */
    static const int numberOfDimensions = p.getNumberOfDimensions();
    double ret = 0; // ret = return value
    for (int q = 0; q < numberOfDimensions; q++)
    {
        ret += p.getPosition()[q] * p.getPosition()[q];
    }
    return ret;
}

double particle_r2(Particle &p1, Particle &p2)
{
    /*
    Calculate (r_1 - r_2)^2 for particle p1 and p2
    */
    static const int numberOfDimensions = p1.getNumberOfDimensions();
    double rdiff;
    double ret = 0;
    for (int q = 0; q < numberOfDimensions; q++)
    {
        rdiff = p1.getPosition()[q] - p2.getPosition()[q];
        ret += rdiff * rdiff;
    }
    return ret;
}

double dot_product(std::vector<double> &v1, std::vector<double> &v2, int numberOfDimensions)
{
    /*
    Computes the dot product of vectors v1 and v2.
    */
    double ret = 0;
    for (int i = 0; i < numberOfDimensions; i++)
        ret += v1[i] * v2[i];

    return ret;
}

void particle_add_rdiff(std::vector<double> &diff, Particle &p1, Particle &p2, double scale)
{
    /*
    Adds a vector term r_1-r_2 to the vector diff. The scale parameter is used since the calculation of "v" has a factor of u'(r_12)/|r_12|
    */
    static const int numberOfDimensions = p1.getNumberOfDimensions();
    for (int i = 0; i < numberOfDimensions; i++)
        diff[i] += scale * (p1.getPosition()[i] - p2.getPosition()[i]);
}

void Particle::saveEquilibrationPosition()
{
    m_EquilibrationPosition = m_position;
}

void Particle::resetEquilibrationPosition()
{
    m_position = m_EquilibrationPosition;
}