#ifndef GLOBAL_H
#define GLOBAL_H

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
	std::vector<double> wall_y_coor;
	std::vector<double> x_coor;
};



#endif
