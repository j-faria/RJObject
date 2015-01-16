#include "MyDistribution.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace DNest3;

/*
We call "parameters" to the hyperparameters of the priors (alpha in the paper)
	perturb_parameters() works on these
	since this leaves the positions unchanged, it will usually not affect the 
	value of the likelihood, and the likelihood will not need to be recomputed

We call "positions" to the parameters of the model, the so-called components x_i
in the paper. 
	

*/



MyDistribution::MyDistribution(double x_min, double x_max,
					double mu_min, double mu_max)
:x_min(x_min)
,x_max(x_max)
,mu_min(mu_min)
,mu_max(mu_max)
{

}

// This is sampling from the prior for the hyperparameters
void MyDistribution::fromPrior()
{
	mu = exp(log(mu_min) + log(mu_max/mu_min)*randomU());
}

// This is a proposal for the hyperparameters
double MyDistribution::perturb_parameters()
{
	double logH = 0.;

	mu = log(mu);
	mu += log(mu_max/mu_min)*pow(10., 1.5 - 6.*randomU())*randn();
	mu = mod(mu - log(mu_min), log(mu_max/mu_min)) + log(mu_min);
	mu = exp(mu);

	return logH;
}

// vec[0] = "position" (log-period)
// vec[1] = amplitude
// vec[2] = phase

double MyDistribution::log_pdf(const std::vector<double>& vec) const
{
	if(vec[0] < x_min || vec[0] > x_max || vec[1] < 0. ||
			vec[2] < 0. || vec[2] > 2.*M_PI)
		return -1E300;

	return -log(mu) - vec[1]/mu;
}

// This is the inverse of the function F, the G function in the paper
void MyDistribution::from_uniform(std::vector<double>& vec) const
{
	vec[0] = x_min + (x_max - x_min)*vec[0];
	vec[1] = -mu*log(1. - vec[1]);
	vec[2] = 2.*M_PI*vec[2];
}

// This is the F (x;α) function from the paper
// it takes a component x and transforms it to a value u that has a uniform 
// distribution between 0 and 1, given α. 
// If the component x consists of a single scalar value, F is the cumulative 
// distribution of the conditional prior
void MyDistribution::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - x_min)/(x_max - x_min);
	vec[1] = 1. - exp(-vec[1]/mu);
	vec[2] = vec[2]/(2.*M_PI);
}

void MyDistribution::print(std::ostream& out) const
{
	out<<mu<<' ';
}

