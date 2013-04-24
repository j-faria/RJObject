#include <ctime>
#include <iostream>
#include <fstream>
#include "RJObject.h"
#include "MassDistributions/Exponential.h"

using namespace DNest3;
using namespace std;

// Demonstration of the MCMC
int main()
{
	// Initialise DNest3's random number generator
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	// Make an object with 2 spatial dimensions, maximum
	// of 100 components, and exponential prior on the masses
	// given a mean, which is log-uniform between 1E-3 and 1E3.
	RJObject<Exponential> r(2, 100, Exponential(1E-3, 1E3));

	// Generate the object from the prior
	r.fromPrior();

	// Open an output file stream to save the samples
	fstream fout("output.txt", ios::out);

	// How many MCMC steps to do
	int steps = 100;

	for(int i=0; i<steps; i++)
	{
		// Make a proposal
		RJObject<Exponential> r2 = r;
		double logH = r2.perturb();

		// Accept the proposal?
		if(randomU() <= exp(logH))
			r = r2;

		r.print(fout);
		fout<<endl;
		cout<<(i+1)<<'/'<<steps<<endl;
	}
	fout.close();
	return 0;
}

