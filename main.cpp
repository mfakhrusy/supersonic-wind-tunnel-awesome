#include "global.hpp"
#include "prandtl_meyer.hpp"

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
	std::vector<double> &k_plus 		= 	parameters.k_plus;
	std::vector<double> &k_minus 		= 	parameters.k_minus;
	std::vector<double> &theta		= 	parameters.theta;
	std::vector<double> &nu			= 	parameters.nu;	//prandtl meyer angle
	std::vector<double> &mach		= 	parameters.mach;
	std::vector<double> &mu			= 	parameters.mu;	//mach angle
	std::vector<double> &area_ratio		= 	parameters.area_ratio;//A/A*
	std::vector<double> &y_coor		=	parameters.y_coor;
	std::vector<double> &wall_y_coor	=	parameters.wall_y_coor;
	std::vector<double> &x_coor		=	parameters.x_coor;

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
	
	//initialize x_coor
	x_coor.push_back(0);
	x_coor.push_back((y_coor[steps] - y_coor[0])/tan(theta[steps]/prandtl_meyer.deg_factor));
	wall_y_coor.push_back(y_coor[0]);
	wall_y_coor.push_back(y_coor[steps]);

	int count = steps;
	//calculate x coordinate
	for (auto i = 0; i < steps; i++) {
		count = count + (steps-i) + 1;
		wall_y_coor.push_back(y_coor[count]);
		if (i == steps-1) {
			x_coor.push_back((y_coor[count] - y_coor[0])/tan(theta[count-1]/prandtl_meyer.deg_factor));
		} else {
			x_coor.push_back((y_coor[count] - y_coor[0])/tan(theta[count]/prandtl_meyer.deg_factor));
		}

	}
	
	//write to file
	prandtl_meyer.write_to_file(parameters, mu.size());
}

