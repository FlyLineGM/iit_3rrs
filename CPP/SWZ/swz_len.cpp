/** @file 3rrs_Main.cpp 
The program is used to find designs of 3rrs which satisfy the kinematic design constraints using the concept of "safe working zone" by searching through the design space of its link lengths.
*/

#include <vector>
#include <complex.h>
#include <iostream>
//#include <time.h>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <fstream>

#include "swz_3rrs.h"
#include "scanner.h"
#include "swz_len.h"

//#define PAR // Comment out if not parallel version
using namespace std;


double motorheight = 0.15; // this is the sum of half of motor height and an offset of 50 mm
double motorwidth = 0.12;
double motorlength = 0.710;


//Main function
//----------------------------------------------------------

// The structure of this program involves looping over a range of the design variables and finding the swz data for each of them. The output is written in radData with all the design variables. 

void swz_len(double des_rb, double des_lcr, double des_lst, double des_rt, double des_del, double XMIN, double XMAX, double YMIN, double YMAX, double X0, double Y0, double RREQ, double Zrange[3])
{
 const int fn_num = 4; // The total number of functions to be iterated for the S1-S4 functions.
 const int param_num = 6; // Number of Design parameters (5+1). z is taken as a parameter to scan each slice

//Initial sample size in the x and y directions set to 0.5 deg and the minimum sample size set to 2 deg.
 const int cSamples = ceil(2*(XMAX-XMIN)*180/M_PI); 

 const int MinSamples = 50;
 
 const double ZRESO = 1.0; // search resolution in Z-direction, in mm.
 const double mTOmm = 1000;// Constant to conver m to mm
 
 const double LARGECONS = 10000;

  double parameters[param_num], radii[fn_num];;
  //double des_x1, des_l1, des_l0, des_r, des_rt, des_lamda, des_del;
  double xmin, ymin, xmax, ymax, zmin, zmax;
  double minRadii;
  int SAMPLESZ=0;

//Scan region: super set of SWZ
  zmax = des_lst+des_lcr;
  zmin = motorheight + sqrt((des_rb+motorwidth)*(des_rb+motorwidth)+ motorlength*motorlength) * tan(15 * M_PI/180);//Min height at which the platform when tilted by 15 deg hits the farthest corner point of the motor. Picture...

  SAMPLESZ = (int)ceil((zmax-zmin)*mTOmm/ZRESO); //numbers of layers to scan
  if(SAMPLESZ<0)
    SAMPLESZ=0;

  //std::cout <<"zvals: "<< zmin << " "<< zmax << " " <<SAMPLESZ<< std::endl;
 
  double printradii[SAMPLESZ+1];


  int SAMPLES[SAMPLESZ+1]; 
  
  int i;
  int i2;

#ifdef PAR
#pragma omp parallel for private(parameters,i,minRadii,radii,xmax,xmin,ymax,ymin, i2) shared(printradii, SAMPLES)//, zmin, zmax, cSamplesX, cSamplesY, SAMPLESZ)
#endif
  for(i=0; i<=SAMPLESZ; i++){
   
    parameters[0] = des_rb; //Radius of base platform
    parameters[1] = des_lcr; //length of crank
    parameters[2] = des_lst; //length of strut
    parameters[3] = des_rt;  // Radius of top platform
    parameters[5] = des_del; // wedge angle

     
    SAMPLES[i] = cSamples;

    xmax = XMAX;
    xmin = XMIN;
    ymax = YMAX;
    ymin = YMIN;
    

//z is the additional parameter for each slice 
    parameters[4]=zmin+(zmax-zmin)*i/SAMPLESZ; //ith z value

    minRadii = LARGECONS;


    for(i2 = 0; i2<fn_num; i2++)
      {
	
      	//cout << X0 << " " << Y0 << " " << xmax << " " << xmin << " " << ymax << " "<< ymin << " " << SAMPLESX[i] << " "<<  SAMPLESY[i]<< " " << MAXREC << " " << minRadii<< endl;
      	
      	radii[i2]=scanner(fn_ptr_selector(i2), X0, Y0, xmax, xmin, ymax, ymin, SAMPLES[i2], parameters);
      	//radii[i2] = scanner(fn_ptr_selector(i2), X0, Y0, XMAX, XMIN, YMAX, YMIN,  SAMPLESX[i], SAMPLESY[i], MAXREC, parameters, i2);

      	if ( minRadii > radii[i2] ) {
      	  minRadii = radii[i2];

      	}
      	
      	printradii[i] = minRadii*180/M_PI;
      	//reducing the search area - heirarchy
      	if(minRadii!=LARGECONS){
      	  xmin = -minRadii;
      	  ymin = -minRadii;
      	  xmax = minRadii;
      	  ymax = minRadii; //??
      	  //ymax = minRadii;
      	  
      	  SAMPLES[i] = ceil(4*minRadii*180/M_PI); //resizing the grid

      	  //prevent grid size becoming too small
      	  if (SAMPLES[i]<MinSamples){
      	    SAMPLES[i] = MinSamples;
      	  }
      	}
      	else {
      	  xmax = XMAX;
      	  xmin = XMIN;
      	  ymax = YMAX;
      	  ymin = YMIN;
      	}
      }
  }//END ONE DESIGN	


  int lencount=0;
  int maxcount = 0;

  Zrange[0] = 0;
  Zrange[1] = zmin;
  Zrange[2] = zmin;

  for(int il = 0; il<SAMPLESZ+1; il++)
  {
    if(printradii[il]>=RREQ)
    {
     lencount++;
    }
    else
    {
      lencount = 0;
    }
    
    if(lencount > maxcount)
    {
      maxcount = lencount;
      Zrange[2] = zmin+(zmax-zmin)*(il)/SAMPLESZ;
      Zrange[1] = zmin+(zmax-zmin)*(il-lencount+1)/SAMPLESZ;
      Zrange[0] = Zrange[2]-Zrange[1];
    } 
  }
}
