#include <iostream>
#include <cmath>
# include "swz_len.h"

int main()
{

  double l[] = {.407,.323,.650,.806};
  double del0 = 0*M_PI/180;
  // SWZ
  double XMIN = -M_PI/6;
  double XMAX = M_PI/6;
  double YMIN = -M_PI/6;
  double YMAX = M_PI/6; 
    
  // center point of the workspace
  double X0 = 0;
  double Y0 = 0;
    
  // radius of the cylinder required in degrees.
  double RREQ = 20;//*M_PI/180;
    
  double zrange[3];

  swz_len(l[0], l[1], l[2], l[3], del0, XMIN, XMAX, YMIN, YMAX, X0, Y0, RREQ, zrange);

  std::cout<<zrange[0]<<std::endl;
  std::cout<<zrange[1]<<std::endl;
  std::cout<<zrange[2]<<std::endl;  
  
	return 0;
}
