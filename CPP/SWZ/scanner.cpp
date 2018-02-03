#include <iostream>	
#include <math.h>

//////////////////////////////////
#define    PI         acos(-1)
#define    EPS        1.0e-10
#define    sign(n)    ((n) < 0) ? -1 : 1
#define    max(A,B)   ((A) > (B)) ? A : B
#define    min(A,B)   ((A) < (B)) ? A : B
//////////////////////////////////
double scanner(double (*function)(double, double, double*), double x0, double y0, double xmax, double xmin, double ymax, double ymin, int samples, double *parameters)
{
	double dmin = min(xmax-xmin, ymax-ymin);
	double r_stepsize = dmin / samples;

	//printf("dmin: %f , r_stepsize: %f\n", dmin,r_stepsize);

	double r;
	double theta;
	double theta_stepsize=PI;
	double s=r_stepsize*theta_stepsize;

	int initialSign;
	double fnValue;
	double x,y;

	x= x0;
	y= y0;

	fnValue=function(x,y,parameters);
	initialSign= sign(fnValue);

//	long fn_samples=0;

	for(r=r_stepsize; r<= dmin/2; r+=r_stepsize)
	{
//If the theta_stepsize is less than 1 deg step_size is not reduced any further
		if(theta_stepsize >= PI/180)
			theta_stepsize=s/r;
//The number of points keeps increasing by the same factor as r.
		for(theta=0; theta < 2*PI; theta+=theta_stepsize)
		{
			x= x0 + r*cos(theta);
			y= y0 + r*sin(theta);

			fnValue=function(x,y,parameters);
			//fn_samples++;
			// std::cout<<x<<'\t'<<y<<'\t'<<fnValue<<std::endl;
			if( (sign(fnValue)) != initialSign || isnan(fnValue) ) {
				return r-r_stepsize;
			}
		}
	}

	return r-r_stepsize;
}
