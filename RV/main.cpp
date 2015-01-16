#include <iostream>
#include "Start.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest3;



int main(int argc, char** argv)
{
	// Load the data
	Data::get_instance().load("fake_data.txt");
	//Data::get_instance().load("14her.rv");
	//return 0;


	MTSampler<MyModel> sampler = setup_mt<MyModel>(argc, argv);
	sampler.run();

	return 0;
}




