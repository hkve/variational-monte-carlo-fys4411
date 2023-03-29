#include <memory>
#include <cmath>
#include <cassert>

#include "interactinggaussian.h"
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"

#include <iostream>
using namespace std;

InteractingGaussian::InteractingGaussian(double alpha, double beta, double a, int num_particles) // a is the parameter in the interaction
{
    assert(alpha >= 0);
    assert(beta == 1 || beta == 2.82843);
    m_interactionTerm = a;               // Parameter in the interaction
    m_numberOfParameters = 2;            // notice for simple gaussian we had 2 also, but now beta will not be 1
    m_numberOfParticles = num_particles; // this is the number of particles in the system. Different from the non interacting case, it is important
                                         // that this is a member variable, since we need to know the number of particles in the system in order to
                                         // compute the interaction term in the quantum force, for example, which do not currently have access to.
    m_parameters.reserve(2);
    m_parameters.push_back(alpha); // m_parameters is the vector of variational parameters
    m_parameters.push_back(beta);
}

double InteractingGaussian::evaluate(std::vector<std::unique_ptr<class Particle>> &particles)
{

    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i]->getPosition()
     * function.
     */
    double r2 = 0; // r2 is the sum of the squared coordinates of the r vector
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);
    double a = m_interactionTerm; // renaming for simplicity in formulas
    double r_ij = 0;
    double interaction = 1;

    for(int i = 0; i < m_numberOfParticles; i++) {
        Particle& particle_i = *particles[i];
        r2 += particle_r2(particle_i);
        r2 += particle_i.getPosition()[2]*particle_i.getPosition()[2]*(beta-1);
    
        for(int j = i+1; j < m_numberOfParticles; j++) {
            Particle& particle_j = *particles[j];
            r_ij = std::sqrt(particle_r2(particle_i, particle_j));
            interaction *= (1 - a/r_ij) * (r_ij > a);
        }
    }
    
    double gaussian = std::exp(-alpha * r2); // Notice this includes beta.
    return gaussian * interaction;
}

double InteractingGaussian::evaluate_w(int proposed_particle_idx, class Particle &proposed_particle, class Particle &old_particle, std::vector<std::unique_ptr<class Particle>> &particles)
{
    /*
     This is the wave function ratio for the Metropolis algorithm.
     It is a clever way to avoid having to evaluate the wave function for all particles at each step.
     The gaussian part is still present, but we also have to recalculate every term where the proposed particle is present (one N product with f_ij)
    */
    static const int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    static const double a = m_interactionTerm;
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);

    double r2_proposed, r2_old;
    r2_proposed = 0;
    r2_old = 0;

    r2_proposed = particle_r2(proposed_particle);
    r2_old = particle_r2(old_particle);

    // beta corrections to r2. Notice this lets us use the same r2, even if beta is not 1
    r2_proposed += (proposed_particle.getPosition()[2] * proposed_particle.getPosition()[2]) * (beta - 1);
    r2_old += (old_particle.getPosition()[2] * old_particle.getPosition()[2]) * (beta - 1);

    double gaussian = std::exp(-2.0 * alpha * (r2_proposed - r2_old)); // Same as non-interactive

    double interaction = 1;
    double r_gj_prime = 0; // |r_g - r_j| in proposed R configuration
    double r_gj = 0;       // |r_g - r_j| in old R configuration
    double delta = 0;      // If any particle distances are less than a, evalute interaction term to 0

    // proposed_idx != i product. Divided into two loops to avoid if statments
    for (int i = 0; i < proposed_particle_idx; i++)
    {
        r_gj_prime = std::sqrt(particle_r2(proposed_particle, *particles[i]));
        r_gj = std::sqrt(particle_r2(old_particle, *particles[i]));
        delta = (r_gj_prime > a) * (r_gj > a);
        if(!delta)
            return 0;
        interaction *= (1.0 - a / r_gj_prime) / (1.0 - a / r_gj); // ratio for relative r_gj distance
    }
    // Same as above but for the indicies after proposed_particle_idx
    for (int i = proposed_particle_idx + 1; i < m_numberOfParticles; i++)
    {
        r_gj_prime = std::sqrt(particle_r2(proposed_particle, *particles[i]));
        r_gj = std::sqrt(particle_r2(old_particle, *particles[i]));
        delta = (r_gj_prime > a) * (r_gj > a);
        if(!delta)
            return 0;
        interaction *= (1.0 - a / r_gj_prime) / (1.0 - a / r_gj);
    }

    return gaussian * interaction * interaction; // Dont forget to square the interaction part :)
}

double InteractingGaussian::computeParamDerivative(std::vector<std::unique_ptr<class Particle>> &particles, int parameterIndex)
{
    /* Note that by derivative, we actually
     * mean the derivative with respect to the variational parameters.
     * The interaction term does not depend on the variational parameters,
     * so we can just return the derivative of the gaussian part wrt alpha.
     * Notice we don't even need to multiply by the interaction term, because we then divide by the WF (TRIPLE CHECK THIS!).
     */
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);

    double r2 = 0;
    double rz2 = 0;
    for (int k = 0; k < m_numberOfParticles; k++)
    {
        Particle &particle = *particles[k];
        r2 += particle_r2(particle);
        // beta correction to r2. Notice this lets us use the same r2, even if beta is not 1
        rz2 = (particle.getPosition()[2] * particle.getPosition()[2]);
        r2 += rz2 * (beta - 1);
    }
    if (parameterIndex == 0)
    {
        return -r2;
    }
    else if (parameterIndex == 1)
    {
        return -alpha * rz2;
    }
    else
    {
        cout << "Error: parameterIndex out of bounds in InteractingGaussian::computeParamDerivative" << endl;
        return 0;
    }
}

double InteractingGaussian::u_p(double r_ij)
{
    /*
    Closed form derivative, u'(r_ij)
    */
    static const double a = m_interactionTerm;
    return a / (r_ij * (r_ij - a));
};

double InteractingGaussian::u_pp(double r_ij)
{
    /*
    Closed form derivative, u''(r_ij)
    */
    static const double a = m_interactionTerm;
    return -(a * (2.0 * r_ij - a)) / (r_ij * r_ij * (r_ij - a) * (r_ij - a));
}

void InteractingGaussian::grad_phi_ratio(std::vector<double> &v, Particle &particle, double alpha, double beta)
{
    /*
    Calculates the ratio of the OB gradient wrt. to "particle"'s coordinate. Divided by OB wf.
    */
    static const int numberOfDimensions = particle.getNumberOfDimensions();
    for (int i = 0; i < numberOfDimensions; i++)
        v[i] = -2 * alpha * particle.getPosition()[i];

    // beta correction to gradient.
    v[2] *= beta;
}

double InteractingGaussian::computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>> &particles)
{
    /* Analitical double derivative. By double derivative, we mean the sum of the Laplacians with respect to the
     * coordinates of each particle. And also we divide by the wave function because it simplifies the
     * calculations in the Metropolis algorithm and it is all we actually need.
     */
    static const int numberOfDimensions = particles.at(0)->getNumberOfDimensions();
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);

    double a = m_interactionTerm; // renaming for simplicity in formulas

    std::vector<double> grad_phi_ratio_k(3); // Storing the OB gradient divided by OB wf
    double r_kj_length;                      // The length of the r_k - r_j vector
    double u_p_kj, u_pp_kj;                  // Store u' and u'' evaluated at r_kj

    double term1 = 0; // First TB term, dot product of v and gradient wf
    double term2 = 0; // Second TB term, simply v dot v
    double term3 = 0; // Third TB term, sum of scalar u' and u''

    double r2_sum_OB = 0; // The sum of all squared coordinates.
    for (int k = 0; k < m_numberOfParticles; k++)
    {
        std::vector<double> v(numberOfDimensions, 0); // For each k, reset v vector

        // first we compute the double derivative of the gaussian part
        Particle &particle_k = *particles[k];

        // Add r_k^2 for OB term
        r2_sum_OB += particle_r2(particle_k);

        // beta correction to r2. Notice this lets us use the same r2, even if beta is not 1
        r2_sum_OB += (particle_k.getPosition()[2] * particle_k.getPosition()[2]) * (beta*beta - 1);

        // Calculate OB gradient and wf ratio
        grad_phi_ratio(grad_phi_ratio_k, particle_k, alpha, beta);

        // Sum over j != k, divded into two parts to avoid if-statement. First do all particles up to k
        for (int j = 0; j < k; j++)
        {
            Particle &particle_j = *particles[j];
            r_kj_length = std::sqrt(particle_r2(particle_k, particle_j)); // r_kj = |r_k - r_j|

            u_p_kj = (r_kj_length > a) ? u_p(r_kj_length) : 0;   // Calculate u'(r_kj) if r_kj is less than a. Else evalute to 0
            u_pp_kj = (r_kj_length > a) ? u_pp(r_kj_length) : 0; // Same for u''(r_kj)

            // Here we add a term of the j != k sum to the vector v. This vector term to add is scaled by a factor u'(r_kj)/r_kj
            // Note that this ADDS to v, not overwriting (cumulative)
            particle_add_rdiff(v, particle_k, particle_j, u_p_kj / r_kj_length);

            // For this j, add the third term constribution
            term3 += u_pp_kj + (2 / r_kj_length) * u_p_kj;
        }
        // The rest of the particles, from k+1 to N. The giblets of the loop is the same.
        for (int j = k + 1; j < m_numberOfParticles; j++)
        {
            Particle &particle_j = *particles[j];
            r_kj_length = std::sqrt(particle_r2(particle_k, particle_j));

            u_p_kj = (r_kj_length > a) ? u_p(r_kj_length) : 0;
            u_pp_kj = (r_kj_length > a) ? u_pp(r_kj_length) : 0;

            particle_add_rdiff(v, particle_k, particle_j, u_p_kj / r_kj_length);

            term3 += u_pp_kj + (2 / r_kj_length) * u_p_kj;
        }

        // Now that we are done with a j != k sum, perform the dot products for first and second term. Add this to the total
        term1 += 2 * dot_product(grad_phi_ratio_k, v, numberOfDimensions);
        term2 += dot_product(v, v, numberOfDimensions);
    }

    // cout << term1 << " " << term2 << " " << term3 << endl;
    double gaussian_double_derivative = 2 * alpha * (2 * alpha * r2_sum_OB - m_numberOfParticles * (numberOfDimensions - 1 + beta)); // The double derivative of the gaussian part
    double interaction_double_derivative = term1 + term2 + term3;                                                                    // Adding the contribution from the three TB terms

    // cout << gaussian_double_derivative << " " << interaction_double_derivative << endl;
    return gaussian_double_derivative + interaction_double_derivative;
}

// Notice that now we need to pass the whole vector of particles, because we need to compute the interaction term.
void InteractingGaussian::quantumForce(std::vector<std::unique_ptr<class Particle>> &particles, Particle &particle, std::vector<double> &force)
{
    static const int numberOfDimensions = particle.getNumberOfDimensions();
    static const double a = m_interactionTerm;
    double alpha = m_parameters.at(0);
    double beta = m_parameters.at(1);

    for (int q = 0; q < numberOfDimensions; q++)
    {
        force.at(q) = -4.0 * alpha * particle.getPosition().at(2); // analytic derivative wrt r_q for the gaussian part
    }

    // beta correction
    force.at(2) += -4.0 * alpha * particle.getPosition().at(2) * (1 - beta);

    // Now we need to add the interaction term
    double r_ij, r_ij_q, norm_rij;
    for (int i = 0; i < m_numberOfParticles; i++)
    {
        Particle &other_particle = *particles.at(i);
        norm_rij = std::sqrt(particle_r2(particle, other_particle));

        if (norm_rij > a) // notice this automatically excludes the particle itself
        {
            for (int q = 0; q < numberOfDimensions; q++)
            {
                r_ij_q = particle.getPosition().at(2) - other_particle.getPosition().at(2);
                force.at(q) += 2 * a * r_ij_q / (norm_rij * norm_rij * (norm_rij - a));
            }
        }
    }
}

void InteractingGaussian::setParameters(std::vector<double> parameters)
{
    assert((int)parameters.size() == m_numberOfParameters);
    m_parameters = parameters;
}