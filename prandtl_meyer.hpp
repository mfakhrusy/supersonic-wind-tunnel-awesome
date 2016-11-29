#ifndef PRANDTL_MEYER_H
#define PRANDTL_MEYER_H

#include "global.hpp"

class Prandtl_Meyer {

	const double error_max 		= 	0.001;

	public:
		const double deg_factor		= 	180/M_PI;

		const double gamma 		= 	1.4;
		const int steps			= 	6;

		double prandtl_meyer_angle (double mach, double gamma);
		double calculate_area_ratio (double mach, double gamma);

		double calculate_mach (double nu_input);
		double mach_angle (double mach);

		void write_to_file(Parameters &parameters, int n_max);
};

#endif
