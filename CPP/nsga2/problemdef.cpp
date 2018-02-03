/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
# include "../SWZ/swz_len.h"
# include "../Dynamics/inertiamodel_3rrs.h"

void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
  double l[4];

  l[0] = xreal[0];//rb  -- Base platform circum radius
  l[1] = xreal[1];//lCr -- Length of crank
  l[2] = xreal[2];//lSt -- Length of strut
  l[3] = xreal[3];//rt  -- Top platform circum radius

  double wAng = xreal[4]*M_PI/180; // Wedge Angle

  // SWZ
  double XMIN = -M_PI/5;
  double XMAX = M_PI/5;
  double YMIN = -M_PI/5;
  double YMAX = M_PI/5; 
    
  // center point of the workspace
  double X0 = 0;
  double Y0 = 0;
    
  // radius of the cylinder required in degrees.
  double RREQ = 20;//*M_PI/180;

  // output data in the order: zrange, zmin, zmax
  double zrange[3];
  swz_len(l[0], l[1], l[2], l[3], wAng, XMIN, XMAX, YMIN, YMAX, X0, Y0, RREQ, zrange);


  // Defining test object
  InertiaModel simulator;
  simulator.setlength(l);

  // Defining payloads
  double mpl[] = {40, 90, 70};//Base of console, human cylinder, Cuboid(Black portion) weight
  double platformdim[] = {1.5, 1.2, 0.02};//l b h of the console base
  double cylrad = 0.25;//Radius of the cylinder
  double cylheight = 1.2;//Height of the cylinder
  double cubside = 0.6;//Side of the cuboid
  double cubheight = 1;//Height of the cuboid
  double cyloffsetX = +0.311;//Position of the cylinder wrt body fixed frame
  double cuboffsetX = -0.114;//Position of the cuboid wrt body fixed frame
  simulator.setpayload(mpl, platformdim, cylrad, cylheight, cyloffsetX, cubside, cubheight, cuboffsetX);


  // scan Zone 1 -- The core region
  double widthZone1 = 0.05;
  double par1[] = { 1, 1}; // Parameters  

  double torqueZone1;
  torqueZone1 = simulator.torquescan(zrange, widthZone1, par1);  

  //scan Zone 2 -- The Outer region
  double widthZone2 = 0.25;

  //Fraction of specified accelerations in the secondary region
  double alpha = xreal[5];
  double beta = xreal[6];

  double par2[] = {alpha, beta};

  double torqueZone2;
  torqueZone2 = simulator.torquescan(zrange, widthZone2, par2);

  // Objectives
  obj[0] = -(alpha + beta); // Maximize sum of parameters

  // Constraints
  double zrangemin = 0.25;
  double maxTorque = 680;
  constr[0] = zrange[0] - zrangemin; // Constraint on minimum z range
  constr[1] = torqueZone1 - maxTorque; 
  constr[2] = torqueZone2 - maxTorque;

  std::cout<<obj[0]<<'\t'<<constr[0]<<'\t'<<constr[1]<<'\t'<<constr[2]<<'\n';
  return;
}
