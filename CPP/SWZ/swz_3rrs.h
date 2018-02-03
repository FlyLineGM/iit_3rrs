/** @file swz_3rrs.h 
The program is used to find the "safe working zone" of a 3rrs-parallel manipulator.
*/
#ifndef SWZ_3rrs_H
#define SWZ_3rrs_H 1
#include<cmath>
//#include "global.h"
#include <iostream>
double alphabeta2c1c2(double alpha, double beta, double &c1, double &c2);
typedef double (*fptr)(double alpha, double beta, double* params);
fptr fn_ptr_selector(int fn_no);

inline void swz_3rrs_rrdyad(double l1, double l2, double x1, double y1, double x2, double y2, double soln[2]);
inline void swz_3rrs_rrdyad1(double l1, double l2, double x1, double y1, double x2, double y2, double &soln);
inline double norm(double vec[3]);
inline double vecAngle(double vec1[3], double vec2[3]);

int swz_3rrs_invkin_theta1(double alpha, double beta, double *params, double &theta1);
int swz_3rrs_invkin_theta2(double alpha, double beta, double *params, double &theta2);
int swz_3rrs_invkin_theta3(double alpha, double beta, double *params, double &theta3);

int swz_3rrs_invkin1(double alpha, double beta, double *params, double sol1[2]);
int swz_3rrs_invkin2(double alpha, double beta, double *params, double sol2[2]);
int swz_3rrs_invkin3(double alpha, double beta, double *params, double sol3[2]);

int swz_3rrs_invkin_delta1(double alpha, double beta, double *params, double &del1);
int swz_3rrs_invkin_delta2(double alpha, double beta, double *params, double &del2);
int swz_3rrs_invkin_delta3(double alpha, double beta, double *params, double &del3);


double swz_3rrs_S11(double alpha, double beta, double* params);
double swz_3rrs_S12(double alpha, double beta, double* params);
double swz_3rrs_S13(double alpha, double beta, double* params);

double swz_3rrs_S2(double alpha, double beta, double* params);

double swz_3rrs_S301(double alpha, double beta, double* params);
double swz_3rrs_S302(double alpha, double beta, double* params);
double swz_3rrs_S303(double alpha, double beta, double* params);
/*
double swz_3rrs_S304(double alpha, double beta, double* params);
double swz_3rrs_S305(double alpha, double beta, double* params);
double swz_3rrs_S306(double alpha, double beta, double* params);

double swz_3rrs_S307(double alpha, double beta, double* params);
double swz_3rrs_S308(double alpha, double beta, double* params);
double swz_3rrs_S309(double alpha, double beta, double* params);
double swz_3rrs_S310(double alpha, double beta, double* params);
double swz_3rrs_S311(double alpha, double beta, double* params);
double swz_3rrs_S312(double alpha, double beta, double* params);

double swz_3rrs_S313(double alpha, double beta, double* params);
double swz_3rrs_S314(double alpha, double beta, double* params);
double swz_3rrs_S315(double alpha, double beta, double* params);
double swz_3rrs_S316(double alpha, double beta, double* params);
double swz_3rrs_S317(double alpha, double beta, double* params);
double swz_3rrs_S318(double alpha, double beta, double* params);
*/

/** swz_3rrs_rrdyad: Inverse Kinematics for 3rrs for finding the angles of the RRDyad of each leg

@param l1 The length of link1 in rrdyad
@param l2 The length of link2 in rrdyad
@param x1 x-position of rrdyad base
@param x2 x-position of rrdyad end
@param y1 y-position of rrdyad base
@param y2 y-position of rrdyad end
@param soln variable to store the solution in the order, theta and gamma

Variables used:
Variable | Description
-------- |------------
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: void
 */
inline void swz_3rrs_rrdyad(double l1, double l2, double x1, double y1, double x2, double y2, double soln[2])
{
  double temp_1, temp_2, temp_3;
  temp_1 = x2-x1;
  temp_2 = y2-y1;
  temp_3 = temp_1*temp_1+temp_2*temp_2;
  soln[0] = -acos((l1*l1+temp_3-l2*l2)/(2*sqrt(temp_3)*l1)) + atan2(temp_2,temp_1);
  soln[1] = atan2(temp_2-l1*sin(soln[0]),temp_1-l1*cos(soln[0]));

}

/** swz_3rrs_rrdyad1: Inverse Kinematics for 3rrs for finding only theta of the RRDyad of each leg

@param l1 The length of link1 in rrdyad
@param l2 The length of link2 in rrdyad
@param x1 x-position of rrdyad base
@param x2 x-position of rrdyad end
@param y1 y-position of rrdyad base
@param y2 y-position of rrdyad end
@param soln variable to store the solution for theta

Variables used:
Variable | Description
@returns: void
 */
inline void swz_3rrs_rrdyad1(double l1, double l2, double x1, double y1, double x2, double y2, double &soln)
{
  double temp_1, temp_2, temp_3;
  temp_1 = x2-x1;
  temp_2 = y2-y1;
  temp_3 = temp_1*temp_1+temp_2*temp_2;
  soln = -acos((l1*l1+temp_3-l2*l2)/(2*sqrt(temp_3)*l1)) + atan2(temp_2,temp_1);
}

inline double norm(double vec[3])
{
  double nVal;
  nVal = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  return nVal;      
}
inline double vecAngle(double vec1[3], double vec2[3])
{
  double ang;
  ang = acos((+vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2])/(norm(vec1)*norm(vec2)));
  return ang;
}
#endif//SWZ_3rrs_H
