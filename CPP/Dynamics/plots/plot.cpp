# include "../inertiamodel_3rrs.h"
# include <fstream>


using namespace std;

void InertiaModel::plot(double zrange[3], double phi[3], double omega[3], double alpha[3] , double zderv[3])
{
	double value;
	double leg = 2;

	std::ofstream myfile;

	// Opening text file to store gravity values
	myfile.open ("gravity_heave_leg3.txt");	


	for(double i=zrange[1]; i<=zrange[2]; i+=0.001)
	{
		zderv[0] = i;

		GqTerm(phi, omega, alpha, zderv, leg, value);

	    // Storing Torque values in text file
	    myfile <<zderv[0]<<'\t'<<value<<'\n';		
	}
    myfile.close();

	// Opening text file to store gravity values
	myfile.open ("mass_heave_leg3.txt");	

	for(double i=zrange[1]; i<=zrange[2]; i+=0.001)
	{
		zderv[0] = i;
		
		MqTerm(phi, omega, alpha, zderv, leg, value);

	    // Storing Torque values in text file
	    myfile <<zderv[0]<<'\t'<<value<<'\n';		
	}
    myfile.close();	

	// Opening text file to store gravity values
	myfile.open ("gravity_phi_leg3.txt");	

	double phiter = 15*M_PI/180;

	zderv[0] = 0.5;

	for(double i=-phiter; i<=phiter; i+=0.005)
	{
		phi[0] = i;

		GqTerm(phi, omega, alpha, zderv, leg, value);

	    // Storing Torque values in text file
	    myfile <<phi[0]<<'\t'<<value<<'\n';		
	}
    myfile.close();

	// Opening text file to store gravity values
	myfile.open ("mass_phi_leg3.txt");	

	for(double i=-phiter; i<=phiter; i+=0.005)
	{
		phi[0] = i;
		
		MqTerm(phi, omega, alpha, zderv, leg, value);

	    // Storing Torque values in text file
	    myfile <<phi[0]<<'\t'<<value<<'\n';		
	}
    myfile.close();	

}