#ifndef _Data_
#define _Data_

#include <vector>

class Data
{
	private:
		std::vector<double> t, y, err;

	public:
		Data();
		void load(const char* filename);

		// Getters
		const std::vector<double>& get_t() const { return t; }
		const std::vector<double>& get_y() const { return y; }
		const std::vector<double>& get_err() const { return err; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

