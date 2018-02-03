// Function to scan the workspace for maximum torque

#include "inertiamodel_3rrs.h"
#include <float.h>

/*
	torquescan function
	Input: z-range, width of the zone, parameters
	Output: Maximum torque in a region of given width (if SWZ range is greater than width)
			Maximum torque inside SWZ (if SWZ range is smaller than width)
			DBL_MAX if SWZ range is zero
*/

double InertiaModel::torquescan(double zrange[3], double width, double par[2])
{
//Find z-location of minimum torque

	double minTorque = DBL_MAX;

	double zmin = 0;
	double temp = -DBL_MAX;

	double salp, calp;


/*
Arrays	       | Inputs
--------------   ----------
Phi            | {Axis 2 variables kx & ky, angle (15 deg)}
omega          | {Axis 2 variables kx & ky, specified velocity (15 deg/sec)}
alpha          | {Axis 2 variables kx & ky, specified acceleration (-250 deg/sec^2)}
zderv          | {position z, specified heave velocity (0.2 m/s), specified heave acceleration (0.8 g)}
*/

/*
Scan

Angular position/velocities/accelerations:: all these are taken about the horizontal axis

--- All poses are covered
--- For velocity only +15 m/s is considered since it contributes less than 1% to torque
--- For Acceleration -250 deg/sec^2

*/
	double phi[3] = {0, 1, 15*M_PI/180};//{kx, ky, tilt angle = 15 deg}
	double omega[3] = {0, 1, 15*M_PI/180};//{kx, ky, tilt speed = 15 deg/sec}
	double alpha[3] = {0, 1, par[0]*-250*M_PI/180};//{kx, ky, tilt acceleration = -250 deg/sec}
	double zderv[3] = {0.5, 0.2, par[1]*0.8*g};//{heave range, heave velocity, tilt speed = 15 deg/sec}

	int ziter = floor(zrange[0]*1000);
	int alpiter = 360;
	int phiiter = 15;

	double torqueVector;


	if(width<zrange[0])
	{

	for(int i = 0; i<ziter; i+=5)
	{
		zderv[0] = zrange[1] + zrange[0]*i/ziter;

		temp = -DBL_MAX;

		// Rotating Tilt axis 360 degrees in steps of 5 degrees
		for(int j=0; j<alpiter; j+= 5)
		{
			calp = cos(j*M_PI/180);
			salp = sin(j*M_PI/180);

			phi[0] = calp; phi[1] = salp;
			omega[0] = calp; omega[1] = salp;
			alpha[0] = calp; alpha[1] = salp;			

			//Tilting between -15 to 15 degrees in steps of 5 deg
			for(int k=-phiiter; k<phiiter; k+=5)
			{
				phi[2] = k*M_PI/180;
				//Compute the max of the 3 actuator torques for the given posn, orntn, vel, accel.
				torquecalcmax(phi, omega, alpha, zderv, torqueVector);

				//Get the maximum torque required for this rotation manuever for a given slice. 
				if(temp < torqueVector)
					temp = torqueVector;
			}
		}

		//If the temp is less than previous min torque, update the mintorque with this value. Since at this slice the toque requirements is less than at the other slices for the same manuever, the z value is stored and the torque is also stored. This is referred as the sweet spot or sweet slice for the manipulator and the core region is chosen to be around this slice.
		if(minTorque > temp)
		{
			minTorque = temp;
			zmin = zderv[0];
		}
	}

	double z[2];

	//Find the maximum torque around that z-location		

//The width is the size of core region that is predefined as a requirement. Trying to fit two half cylinders above and below the zmin plane, if that is within the allowable z range. If not, then cylinders are taken from one of the extremes as per the condition.

	if(zmin - width/2 < zrange[1])
	{
	z[0] = zrange[1] ;
	z[1] = zrange[1] + width;
	}
	else if(zmin + width/2 > zrange[2])
	{
	z[1] = zrange[2];
	z[0] = z[1] - width;
	}
	else
	{
	z[0] = zmin - width/2;
	z[1] = z[0] + width;
	}	

	// std::cout<<z[0]<<'\t'<<z[1]<<'\t';
	// std::cin.ignore();

	ziter = floor((z[1]-z[0])*1000);//The iterative steps is taken to be as many as zrange*1000

	temp = -DBL_MAX;

//Maximum torque is computed inside the core region 
//	#pragma omp parallel for private(phi, omega, alpha, salp, calp, zderv) shared(temp)
	for(int i = 0; i<ziter; i+=5)
	{
		zderv[0] = z[0] + (z[1]-z[0]) * i /ziter;

		for(int j=0; j<alpiter; j+= 5)
		{
			calp = cos(j*M_PI/180);
			salp = sin(j*M_PI/180);

		// Rotating Tilt axis 360 degrees in steps of 5 degrees
			phi[0] = calp; phi[1] = salp;
			omega[0] = calp; omega[1] = salp;
			alpha[0] = calp; alpha[1] = salp;			

			for(double k=-phiiter; k<phiiter; k+=2.5)
			{
				phi[2] = k*M_PI/180;
				torquecalcmax(phi, omega, alpha, zderv, torqueVector);
				if(temp < torqueVector)
					temp = torqueVector;
			}	

		}
	}
	return temp;
	}
	else if(zrange[0]==0)//If the range of z is 0
		return DBL_MAX;
	else//Case when the zrange and the core range happen to be equal:: then the scanning is done for the entire region and max torque is returned.
	{

	ziter = floor(zrange[0]*1000);

	temp = -DBL_MAX;

	for(int i = 0; i<ziter; i+=5)
	{
		zderv[0] = zrange[1] + zrange[0] * i /ziter;

		// Rotating Tilt axis 360 degrees in steps of 5 degrees
		for(int j=0; j<alpiter; j+= 5)
		{
			calp = cos(j*M_PI/180);
			salp = sin(j*M_PI/180);

			phi[0] = calp; phi[1] = salp;
			omega[0] = calp; omega[1] = salp;
			alpha[0] = calp; alpha[1] = salp;			

			for(int k=-phiiter; k<phiiter; k+=5)
			{
				phi[2] = k*M_PI/180;
				torquecalcmax(phi, omega, alpha, zderv, torqueVector);
				if(temp < torqueVector)
					temp = torqueVector;
			}	

		}
	}
	return temp;
	}
}
