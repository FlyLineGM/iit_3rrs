# include "../inertiamodel_3rrs.h"
# include "../../SWZ/swz_len.h"
# include <float.h>

using namespace std;

int main()
{

  double l[] = {0.300, 0.284, 0.535, 0.518};
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

  double zrange[] = {0.249, 0.372, 0.621 };

  swz_len(l[0], l[1], l[2], l[3], del0, XMIN, XMAX, YMIN, YMAX, X0, Y0, RREQ, zrange);

  cout<<zrange[0]<<'\t'<<zrange[1]<<'\t'<<zrange[2]<<'\n';
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

  //Plotting
  double phi[3] = {0, 1, 15*M_PI/180};
  double omega[3] = {0, 1, 15*M_PI/180};
  double alpha[3] = {0, 1, -250*M_PI/180};
  double zderv[3] = {0.5, 0.2, 0.8*9.81};

  test.plot(zrange, phi, omega, alpha, zderv);

 }  