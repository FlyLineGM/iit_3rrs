# include "inertiamodel_3rrs.h"
# include <float.h>
# include "../SWZ/swz_len.h"

using namespace std;

int main()
{

  double l[] = {.407,.323,.650,.806};
  double del0 = 0* M_PI/180; 

  // SWZ
  // search space for alpha and beta
  double XMIN = -M_PI/5;
  double XMAX = M_PI/5;
  double YMIN = -M_PI/5;
  double YMAX = M_PI/5; 
    
  // center point of the workspace
  double X0 = 0;
  double Y0 = 0;
    
  // radius of the cylinder required in degrees.
  double RREQ = 20;//*M_PI/180;

  double zrange[] = {0.274, 0.418, 0.693 };

  swz_len(l[0], l[1], l[2], l[3], del0, XMIN, XMAX, YMIN, YMAX, X0, Y0, RREQ, zrange);
 
  // cout<<zrange[2]<<endl;
  // Defining test object
	InertiaModel test;
  test.setlength(l);

  // Defining payloads
  double mpl[] = {40, 90, 70};
  double platformdim[] = {1.5, 1.2, 0.02};
  double cylrad = 0.25;
  double cylheight = 1.2;
  double cubside = 0.6;
  double cubheight = 1;
  double cyloffsetX = 0.20;
  double cuboffsetX = -0.55;
  test.setpayload(mpl, platformdim, cylrad, cylheight, cyloffsetX, cubside, cubheight, cuboffsetX);

/*
  // Inverse Dynamics
  double phi[3] = {0, 1, 8*M_PI/180};
  double omega[3] = {0, 1, 15*M_PI/180};
  double alpha[3] = {0, 1, -250*M_PI/180};
  double zderv[3] = {0.61, 0.2, 0.8*9.81};

  double torque;
  test.torquecalcmax(phi, omega, alpha, zderv, torque);

  cout<<torque<<endl;
*/

  //Dynamic indices scan

 
//  double indices[2];
//  test.dynScan(zrange,indices);
//  std::cout<<indices[0]<<'\t'<<indices[1]<<std::endl;

/*  
  VectorXd q(9);
  double eig[3];

  q<< -0.2763, 0.0544, -0.7818, 1.5059, 1.5576, 1.4224, 0.1135, -0.0099, 0.5486;

  test.calcMthetaEigen(q, eig);

  std::cout<<eig[0]<<'\t'<<eig[1]<<'\t'<<eig[2]<<'\t';
*/

// scan the zone
  double width = 0.20;
  double par[] = { 1, 1}; // Parameters  

  double Torque;
  Torque = test.torquescan(zrange, width, par);
  cout<<Torque<<endl;

  for(int i=0; i<3; i++)
  {
    cout<<zrange[i]<<endl;
  }

  // test.output();
 /*
  //Plotting
  double phi[3] = {0, 1, 15*M_PI/180};
  double omega[3] = {0, 1, 15*M_PI/180};
  double alpha[3] = {0, 1, -250*M_PI/180};
  double zderv[3] = {0.5, 0.2, 0.8*9.81};

  test.plot(zrange, phi, omega, alpha, zderv);
*/
  return 0;
}