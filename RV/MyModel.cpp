#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include "radial_v.h"
#include <cmath>
#include <ctime>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
/* num_dimensions, max_num_components, fixed, dist */
:objects(3, 1, false, MyDistribution(-10., 10., 1E-3, 1E3)), mu(Data::get_instance().get_t().size())
{

}

void MyModel::fromPrior()
{
	objects.fromPrior();
	objects.consolidate_diff();
	sigma = exp(log(1E-3) + log(1E6)*randomU());
	calculate_mu();
}

void MyModel::calculate_mu()
{
	// Get the times from the data
	const vector<double>& t = Data::get_instance().get_t();

	
	std::vector<double> params;

	// Update or from scratch?
	bool update = (objects.get_added().size() < objects.get_components().size());
	//bool update = true;

	// Get the components
	const vector< vector<double> >& components = (update)?(objects.get_added()):
				(objects.get_components());
	//const vector< vector<double> >& components = objects.get_components();

	
	//if(components.size() != 0)
		//std::cout << components.size() << " " << components[0].size() << "\n";

	// Zero the signal
	if(!update)
		mu.assign(mu.size(), 0.);

	double prd,amp,ecc,omega,chi,vsys;
	for(size_t j=0; j<components.size(); j++) // cycle through the components (the planets)
	{
		prd = exp(components[j][0]);
		amp = components[j][1];
		//ecc = components[j][2];
		ecc = 0.3;
		//omega = components[j][3];
		omega = components[j][2];
		//chi = components[j][4];
		chi = 0.;
		//vsys = components[j][5];
		vsys = 0.;
		//std::cout << prd << amp << ecc << omega << chi << vsys << "\n";
		for(size_t i=0; i<t.size(); i++) // cycle through the data points
			mu[i] += Suto::rad_v(t[i], prd,amp,ecc,omega,chi,vsys, 1 );
	}

/*	double T, A, phi;
	for(size_t j=0; j<components.size(); j++) // cycle through the components (the planets)
	{
		T = exp(components[j][0]);
		A = components[j][1];
		phi = components[j][2];
		for(size_t i=0; i<t.size(); i++) // cycle through the data points
			mu[i] += A*sin(2.*M_PI*t[i]/T + phi);
	}
*/}

double MyModel::perturb()
{
	double logH = 0.;

	if(randomU() <= 0.75)
	{
		logH += objects.perturb();
		objects.consolidate_diff();
		calculate_mu();
	}
	else
	{
		sigma = log(sigma);
		sigma += log(1E6)*randh();
		sigma = mod(sigma - log(1E-3), log(1E6)) + log(1E-3);
		sigma = exp(sigma);
	}

	return logH;
}

double MyModel::logLikelihood() const
{
	// Get the data
	const vector<double>& y = Data::get_instance().get_y();

	double logL = 0.;
	double var = sigma*sigma;
	for(size_t i=0; i<y.size(); i++)
		logL += -0.5*log(2.*M_PI*var) - 0.5*pow(y[i] - mu[i], 2)/var;

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
	out<<sigma<<' ';
	objects.print(out); out<<' ';
}

string MyModel::description() const
{
	return string("objects, sigma");
}

