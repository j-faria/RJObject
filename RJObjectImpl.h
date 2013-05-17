#include <cmath>
#include <iostream>
#include <iomanip>
#include "RandomNumberGenerator.h"
#include "Utils.h"

template<class SpatialDist, class MassDist>
RJObject<SpatialDist, MassDist>::RJObject(int num_dimensions, int max_num_components,
				const SpatialDist& spatial_dist,
				const MassDist& mass_dist)
:num_dimensions(num_dimensions)
,max_num_components(max_num_components)
,spatial_dist(spatial_dist)
,mass_dist(mass_dist)
{

}

template<class SpatialDist, class MassDist>
void RJObject<SpatialDist, MassDist>::fromPrior()
{
	// Generate the spatial distribution parameters
	spatial_dist.fromPrior();

	// Generate the mass distribution parameters
	mass_dist.fromPrior();

	// Generate from {0, 1, 2, ..., max_num_components}
	num_components = DNest3::randInt(max_num_components + 1);

	// Resize the vectors of positions and masses
	masses.resize(num_components);
	u_masses.resize(num_components);

	u_positions.resize(num_components, std::vector<double>(num_dimensions));
	positions.resize(num_components, std::vector<double>(num_dimensions));

	// Generate positions and masses
	for(int i=0; i<num_components; i++)
	{
		for(int j=0; j<num_dimensions; j++)
		{
			u_positions[i][j] = DNest3::randomU();
			positions[i][j] = u_positions[i][j];
		}
		spatial_dist.position_from_uniform(positions[i]);

		u_masses[i] = DNest3::randomU();
		masses[i] = mass_dist.mass_cdf_inv(u_masses[i]);
	}
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::perturb_masses(double chance, double scale)
{
	if(num_components == 0)
		return 0.;

	// A flag for whether each component gets changed or not
	std::vector<bool> change(num_components, false);
	int count = 0;
	for(int i=0; i<num_components; i++)
	{
		if(DNest3::randomU() <= chance)
		{
			change[i] = true;
			count++;
		}
	}
	// At least do one...
	if(count == 0)
		change[DNest3::randInt(num_components)] = true;

	for(int i=0; i<num_components; i++)
	{
		if(change[i])
		{
			// Perturb (in uniform coordinate system)
			u_masses[i] += scale*DNest3::randn();
			u_masses[i] = DNest3::mod(u_masses[i], 1.);

			// Transform
			masses[i] = mass_dist.mass_cdf_inv(u_masses[i]);
		}
	}

	return 0.;
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::perturb_positions(double chance, double scale)
{
	if(num_components == 0)
		return 0.;

	// A flag for whether each component gets changed or not
	std::vector<bool> change(num_components, false);
	int count = 0;
	for(int i=0; i<num_components; i++)
	{
		if(DNest3::randomU() <= chance)
		{
			change[i] = true;
			count++;
		}
	}
	// At least do one...
	if(count == 0)
		change[DNest3::randInt(num_components)] = true;

	for(int i=0; i<num_components; i++)
	{
		if(change[i])
		{
			// Perturb
			for(int j=0; j<num_dimensions; j++)
			{
				u_positions[i][j] += scale*DNest3::randn();
				u_positions[i][j] = DNest3::mod(u_positions[i][j],
								1.);
				positions[i][j] = u_positions[i][j];
			}
			spatial_dist.position_from_uniform(positions[i]);
		}
	}

	return 0.;
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::add_component()
{
	if(num_components >= max_num_components)
	{
		std::cerr<<"# WARNING: Trying to add component ";
		std::cerr<<"but already full!"<<std::endl;
		return 0.;
	}

	// Increment counter
	num_components++;

	// Generate position
	std::vector<double> pos(num_dimensions);
	for(int j=0; j<num_dimensions; j++)
		pos[j] = DNest3::randomU();
	u_positions.push_back(pos);
	spatial_dist.position_from_uniform(pos);
	positions.push_back(pos);

	// Generate mass
	u_masses.push_back(DNest3::randomU());
	masses.push_back(mass_dist.mass_cdf_inv(u_masses.back()));

	return 0.;
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::perturb_num_components(double scale)
{
	double logH = 0.;

	// Work out how many components we will have after the change
	double delta = max_num_components*scale*DNest3::randn();
	int diff = (int)floor(delta);
	// In case diff is zero, make it + 1
	if(diff == 0)
	{
		if(DNest3::randomU() <= 0.5)
			diff = -1;
		else
			diff =  1;
	}
	int new_num_components = DNest3::mod(num_components + diff, max_num_components + 1);
	diff = new_num_components - num_components;

	// Now do the required changes
	if(diff > 0)
	{
		for(int i=0; i<diff; i++)
			logH += add_component();
	}
	else
	{
		for(int i=0; i<-diff; i++)
			logH += remove_component();
	}

	return logH;
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::perturb()
{
	double logH = 0.;

	int which = DNest3::randInt(5);

	if(which == 0)
	{
		// Do some birth or death
		logH += perturb_num_components(
				pow(10., 1.5 - 6.*DNest3::randomU()));
	}
	else if(which == 1)
	{
		// Change the spatial distribution parameters
		if(DNest3::randomU() <= 0.5)
			logH += spatial_dist.perturb1(u_positions, positions);
		else
			logH += spatial_dist.perturb2(u_positions, positions);
	}
	else if(which == 2)
	{
		logH += perturb_positions(pow(10., 0.5 - 4.*DNest3::randomU()),
					pow(10., 1.5 - 6.*DNest3::randomU()));
	}
	else if(which == 3)
	{
		// Change the mass distribution parameters
		if(DNest3::randomU() <= 0.5)
			logH += mass_dist.perturb1(u_masses, masses);
		else
			logH += mass_dist.perturb2(u_masses, masses);
	}
	else if(which == 4)
	{
		logH += perturb_masses(pow(10., 0.5 - 4.*DNest3::randomU()),
					pow(10., 1.5 - 6.*DNest3::randomU()));
	}

	return logH;
}

template<class SpatialDist, class MassDist>
double RJObject<SpatialDist, MassDist>::remove_component()
{
	if(num_components <= 0)
	{
		std::cerr<<"# WARNING: Trying to remove component ";
		std::cerr<<"but already empty!"<<std::endl;
		return 0.;
	}

	// Find one to delete
	int i = DNest3::randInt(num_components);

	// Delete mass
	u_masses.erase(u_masses.begin() + i);
	masses.erase(masses.begin() + i);

	// Delete position
	u_positions.erase(u_positions.begin() + i);
	positions.erase(positions.begin() + i);

	// Decrement counter
	num_components--;

	return 0.;
}

template<class SpatialDist, class MassDist>
void RJObject<SpatialDist, MassDist>::print(std::ostream& out)
{
	out<<std::setprecision(12);
	out<<num_dimensions<<' '<<max_num_components<<' ';
	spatial_dist.print(out); out<<' ';
	mass_dist.print(out); out<<' ';
	out<<num_components<<' ';

	// Write out masses
	for(int i=0; i<num_components; i++)
		out<<masses[i]<<' ';

	// Pad with zeros (turned-off components)
	for(int i=num_components; i<max_num_components; i++)
		out<<0.<<' ';

	// Write out positions (all of first coordinate,
	// then all of second coordinate, etc)
	for(int j=0; j<num_dimensions; j++)
	{
		for(int i=0; i<num_components; i++)
			out<<positions[i][j]<<' ';

		// Pad with zeros (turned-off components)
		for(int i=num_components; i<max_num_components; i++)
			out<<0.<<' ';
	}
}

