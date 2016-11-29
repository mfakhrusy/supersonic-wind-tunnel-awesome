#include "global.hpp"
#include "prandtl_meyer.hpp"

double Prandtl_Meyer::prandtl_meyer_angle(double mach, double gamma) {
	
	double a = sqrt((gamma + 1)/(gamma - 1));
	double b = atan(sqrt((gamma - 1)/(gamma + 1)*(pow(mach,2) - 1)))*deg_factor;
	double c = atan(sqrt(pow(mach,2) - 1))*deg_factor;

	return a*b-c;
}

double Prandtl_Meyer::calculate_mach(double nu_input) {

	double mach_comp = 1.0;
	double error;

	do {
		double nu_comp 	= prandtl_meyer_angle(mach_comp, gamma);	
		error = fabs(nu_input - nu_comp);
		mach_comp = mach_comp + 0.00001;

	} while(error > error_max);

	return mach_comp;
}

double Prandtl_Meyer::mach_angle(double mach) {
	return asin(1/mach)*deg_factor;
}

double Prandtl_Meyer::calculate_area_ratio (double mach, 
		double gamma) {
	double the_power 	= (gamma + 1)/(2*(gamma - 1));
	double a 		= 2/(gamma + 1);
	double b		= 1 + pow(mach,2)*(gamma - 1)*0.5;

	return (1/mach)*pow((a*b),the_power);
}

void Prandtl_Meyer::write_to_file(Parameters &parameters, 
		int n_max) {

	std::ofstream output_file;
	output_file.open ("output/Parameters_Output.dat");
	for (auto i = 0; i < n_max; i++) {
		output_file << i << " ";
		output_file << parameters.k_minus[i] << " ";
		output_file << parameters.k_plus[i] << " ";
		output_file << parameters.theta[i] << " ";
		output_file << parameters.nu[i] << " ";
		output_file << parameters.mach[i] << " ";
		output_file << parameters.mu[i] << " ";
		output_file << parameters.area_ratio[i] << " ";
		output_file << parameters.y_coor[i] << " ";
		output_file << std::endl;
	}
	output_file.close();

	output_file.open("output/coordinate.dat");
	for (auto i = 0; i <= steps+1; i++) {
		output_file << parameters.x_coor[i] << " ";
		output_file << parameters.wall_y_coor[i] << " ";
		output_file << std::endl;
	}
	output_file.close();
}


