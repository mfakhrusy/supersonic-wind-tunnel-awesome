#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

struct Parameters {
	
	//declare + define vector of computation
	std::vector<double> k_plus;
	std::vector<double> k_minus;
	std::vector<double> theta;
	std::vector<double> nu;			//prandtl meyer angle
	std::vector<double> mach;
	std::vector<double> mu;			//mach angle
	std::vector<double> area_ratio;
	std::vector<double> y_coor;
	std::vector<double> x_coor;
};

class Prandtl_Meyer {

	const double error_max 		= 	0.001;
	double deg_factor 		= 	180/M_PI;

	public:

		const double gamma 		= 	1.4;
		const int steps			= 	10;

		double prandtl_meyer_angle (double mach, double gamma);
		double calculate_area_ratio (double mach, double gamma);

		double calculate_mach (double nu_input);
		double mach_angle (double mach);

		void write_to_file(Parameters &parameters, int n_max);
};

int main() {

	Prandtl_Meyer prandtl_meyer;
	Parameters parameters;

	//declare parameters
	const double gamma 		=	prandtl_meyer.gamma;
	
	//declare input
	const int steps			= 	prandtl_meyer.steps;
	const double M_final		= 	1.8;
	const double theta_init 	= 	0.375;
	const double nu_init 		= 	0.375;
	const double test_area		=	0.0169; 	//m^2

	//calculate maximum theta
	const double theta_max 		= 	prandtl_meyer.prandtl_meyer_angle(M_final, gamma)*0.5;

	//calculate theta steps for each step
	const double theta_steps 	=	(theta_max-theta_init)/steps;

	//calculate throat area
	const double throat_area	=	test_area/prandtl_meyer.calculate_area_ratio(M_final, gamma);

	//declare + define vector of computation
	std::vector<double> &k_plus 	= 	parameters.k_plus;
	std::vector<double> &k_minus 	= 	parameters.k_minus;
	std::vector<double> &theta	= 	parameters.theta;
	std::vector<double> &nu		= 	parameters.nu;	//prandtl meyer angle
	std::vector<double> &mach	= 	parameters.mach;
	std::vector<double> &mu		= 	parameters.mu;	//mach angle
	std::vector<double> &area_ratio	= 	parameters.area_ratio;//A/A*
	std::vector<double> &y_coor	=	parameters.y_coor;
	std::vector<double> &x_coor	=	parameters.x_coor;

	//initialize both k, theta and nu value
	theta.push_back(theta_init);
	nu.push_back(nu_init);
	k_minus.push_back(theta[0] + nu[0]);
	k_plus.push_back(theta[0] - nu[0]);
	for (auto i = 1; i <= steps; i++) {
		theta.push_back(theta_init + i*theta_steps);
		nu.push_back(nu_init + i*theta_steps);	

		k_minus.push_back(theta[i] + nu[i]);
		k_plus.push_back(theta[i] - nu[i]);
	}

	//give value to steps + 1 index
	theta.push_back(theta.back());
	nu.push_back(nu.back());
	k_minus.push_back(k_minus.back());
	k_plus.push_back(k_plus.back());

	//initialize x_coor
	x_coor.push_back(0);

	//start looping!~~~~!~~~!
	for (auto n = steps; n > 0; n--) {
		
		for (auto i = 0; i < n; i++) {
			theta.push_back(0 + i*theta_steps);
			nu.push_back(k_minus[steps - n + 1] + i*theta_steps);	

			k_minus.push_back(theta.back() + nu.back());
			k_plus.push_back(theta.back() - nu.back());
		}

		theta.push_back(theta.back());
		nu.push_back(nu.back());
		k_minus.push_back(k_minus.back());
		k_plus.push_back(k_plus.back());

	}


	//looping for calculating mach number and mach angle, area ratio, y_coor
	for (auto i = 0; i<k_minus.size(); i++) {
		mach.push_back(prandtl_meyer.calculate_mach(nu[i]));
		mu.push_back(prandtl_meyer.mach_angle(mach[i]));
		area_ratio.push_back(prandtl_meyer.calculate_area_ratio(mach[i],gamma));
		y_coor.push_back(0.5*sqrt(throat_area*area_ratio[i]));
	}
	
	for (auto i = 0; i <= steps; i++) {
	}
	
	//write to file
	prandtl_meyer.write_to_file(parameters, mu.size());

}

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
	output_file.open ("Parameters_Output.dat");
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
	
}

