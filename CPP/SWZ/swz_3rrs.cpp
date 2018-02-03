/* Notes:

Update- April 19, 2015

The S3 functions have been changed. Previously they just gave the boundary at which the link interference changes. Correspondingly scanner has also been modified to incorporate exact "0"s for boolean functions.


@file swz_3rrs.cpp 
The program is used to find the "safe working zone" of a 3rrs-parallel manipulator.
*/

//#include "global.h"
#include "swz_3rrs.h"
#include <math.h>

const double SQRT3 = sqrt(3);
double DELTA_LIMIT_MAX = 40*M_PI/180;

using namespace std;
/** fn_ptr_selector: Function to select a desired function for the SWZ of 3rrs

@param fn_no Number of the desired function

Variables used:
Variable | Description
-------- |------------
fn       | An array to hold the pointers to the functions

@returns: pointer to the desired function
 */
fptr fn_ptr_selector(int fn_no)
{

  double (*fn[7])(double,double,double*);

  fn[0] = swz_3rrs_S11;
  fn[1] = swz_3rrs_S12;
  fn[2] = swz_3rrs_S13;
  fn[3] = swz_3rrs_S2;

  fn[4] = swz_3rrs_S301;
  fn[5] = swz_3rrs_S302;
  fn[6] = swz_3rrs_S303;
/*
  fn[10] = swz_3rrs_S307;
  fn[11] = swz_3rrs_S308;
  fn[12] = swz_3rrs_S309;
  fn[13] = swz_3rrs_S310;
  fn[14] = swz_3rrs_S311;
  fn[15] = swz_3rrs_S312;
  */

  return fn[fn_no];
}



/** swz_3rrs_invkin_theta1: Inverse Kinematics to find theta1, crank angle of 1st leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param theta1 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_theta1(double alpha, double beta, double *params, double &theta1)
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
      
  double calpsqr = c_alpha*c_alpha;
  double salpsqr = s_alpha * s_alpha;
  double cbetasqr =  c_beta * c_beta;
  double sbetasqr = s_beta * s_beta;  

  double q1[3]={-((ra*(calpsqr - 2*c_alpha*c_beta - 3*cbetasqr + salpsqr*sbetasqr))/(2 + 2*c_alpha*c_beta)),0,z - ra*s_beta};
  swz_3rrs_rrdyad1(lcr, lst, rb, 0, q1[0], q1[2], theta1);

  return 0;
}

/** swz_3rrs_invkin_theta2: Inverse Kinematics to find theta1, crank angle of 2nd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param theta2 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_theta2(double alpha, double beta, double *params, double &theta2)
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
  
  double q2[3]={(3*ra)/(4 + 4*c_alpha*c_beta) + (ra*cos(2*alpha))/(4 + \
4*c_alpha*c_beta) + (4*ra*c_alpha*c_beta)/(4 + 4*c_alpha*c_beta) - \
(2*ra*cos(2*beta)*s_alpha*s_alpha)/(4 + 4*c_alpha*c_beta) - \
(2*SQRT3*ra*s_alpha*sin(2*beta))/(4 + 4*c_alpha*c_beta),0,(2*z + SQRT3*ra*c_beta*s_alpha + ra*s_beta)/2.};
  
  swz_3rrs_rrdyad1(lcr, lst, rb, 0, q2[0], q2[2], theta2);

  return 0;
}

/** swz_3rrs_invkin_theta3: Inverse Kinematics to find theta3, crank angle of 3rd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param theta3 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_theta3(double alpha, double beta, double *params, double &theta3)
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
  
  
  double q3[3]={(ra*(3 + cos(2*alpha) + 4*c_alpha*c_beta - \
2*cos(2*beta)*s_alpha*s_alpha + 2*SQRT3*s_alpha*sin(2*beta)))/(4 + \
4*c_alpha*c_beta),0,z - (SQRT3*ra*c_beta*s_alpha)/2. + (ra*s_beta)/2.};
  
  swz_3rrs_rrdyad1(lcr, lst, rb, 0, q3[0], q3[2], theta3);

  return 0;
}

/** swz_3rrs_invkin_gamma1: Inverse Kinematics to find gamma1, strut angle of 1st leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param gamma1 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly
soln     | Solution to store the inverse kinematics solution

@returns: 0, if solved properly

*/
int swz_3rrs_invkin1(double alpha, double beta, double *params, double sol1[2])
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
    
  double calpsqr = c_alpha*c_alpha;
  double salpsqr = s_alpha * s_alpha;
  double cbetasqr =  c_beta * c_beta;
  double sbetasqr = s_beta * s_beta;

  double q1[3]={-((ra*(calpsqr - 2*c_alpha*c_beta - 3*cbetasqr + salpsqr*sbetasqr))/(2 + 2*c_alpha*c_beta)),0,z - ra*s_beta};
  swz_3rrs_rrdyad(lcr, lst, rb, 0, q1[0], q1[2], sol1);

  return 0;
}

/** swz_3rrs_invkin_gamma2: Inverse Kinematics to find gamma2, strut angle of 2nd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param gamma2 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly
soln     | Solution to store the inverse kinematics solution

@returns: 0, if solved properly

*/
int swz_3rrs_invkin2(double alpha, double beta, double *params, double sol2[2])
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
  
  
  double q2[3]={(3*ra)/(4 + 4*c_alpha*c_beta) + (ra*cos(2*alpha))/(4 + \
4*c_alpha*c_beta) + (4*ra*c_alpha*c_beta)/(4 + 4*c_alpha*c_beta) - \
(2*ra*cos(2*beta)*s_alpha*s_alpha)/(4 + 4*c_alpha*c_beta) - \
(2*SQRT3*ra*s_alpha*sin(2*beta))/(4 + 4*c_alpha*c_beta),0,(2*z + SQRT3*ra*c_beta*s_alpha + ra*s_beta)/2.};

  
  swz_3rrs_rrdyad(lcr, lst, rb, 0, q2[0], q2[2], sol2);

  return 0;
}

/** swz_3rrs_invkin_gamma3: Inverse Kinematics to find gamma3, strut angle of 3rd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param gamma3 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b        | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly
soln     | Solution to store the inverse kinematics solution

@returns: 0, if solved properly

*/
int swz_3rrs_invkin3(double alpha, double beta, double *params, double sol3[2])
{
  double rb,lcr,lst,ra,z;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  z = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_beta = cos(beta), s_beta = sin(beta);
  
  
  double q3[3]={(ra*(3 + cos(2*alpha) + 4*c_alpha*c_beta - \
2*cos(2*beta)*s_alpha*s_alpha + 2*SQRT3*s_alpha*sin(2*beta)))/(4 + \
4*c_alpha*c_beta),0,z - (SQRT3*ra*c_beta*s_alpha)/2. + (ra*s_beta)/2.};
  
  swz_3rrs_rrdyad(lcr, lst, rb, 0, q3[0], q3[2], sol3);

  return 0;
}

/** swz_3rrs_invkin_delta1: Inverse Kinematics to find delta1, spherical joint angle of 1st leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param delta1 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_delta1(double alpha, double beta, double *params, double &del1)
{
  double del0;
  double sol1[2], gamma1;
  double vecWedge1[3], stVec1[3];

  del0 = params[5];

  double s_alpha = sin(alpha);
  double c_alpha = cos(alpha);
  double s_beta = sin(beta);
  double c_beta = cos(beta);

  double sdel = sin(del0);
  double cdel = cos(del0);

  double salpsqr = s_alpha * s_alpha;

  if(swz_3rrs_invkin1(alpha, beta, params, sol1)) 
    std::cout<<"Error: Check Inverse Kinematics";

  gamma1 = sol1[1];

  // Direction cosines of strut vector
  stVec1[0] = cos(gamma1);
  stVec1[1] = 0;
  stVec1[2] = sin(gamma1);

  double temp1 = 1 + c_alpha*c_beta;
  double temp = sqrt(temp1 * temp1);

  // Direction cosines of wedge vector
  vecWedge1[0] = cdel*s_beta - (c_beta*(c_alpha + c_beta)*sdel)/
    temp;
  vecWedge1[1] = -(c_beta*cdel*s_alpha) - (-((c_alpha*s_alpha*s_beta)/
         temp) + 
      ((c_alpha + c_beta)*s_alpha*s_beta)/temp)*
    sdel;
  vecWedge1[2] = c_alpha*c_beta*cdel - (-((c_alpha*(c_alpha + c_beta)*s_beta)/
         temp) - 
      (salpsqr*s_beta)/temp)*sdel;   

  del1 = vecAngle(vecWedge1, stVec1);
//  std::cout<<del1<<std::endl;
  return 0;
}

/** swz_3rrs_invkin_delta2: Inverse Kinematics to find delta1, spherical joint angle of 2nd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param delta2 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b      | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_delta2(double alpha, double beta, double *params, double &del2)
{
  double del0;
  double sol2[2], gamma2;
  double vecWedge2[3], stVec2[3];

  del0 = params[5];

  double s_alpha = sin(alpha);
  double c_alpha = cos(alpha);
  double s_beta = sin(beta);
  double c_beta = cos(beta);

  double sdel = sin(del0);
  double cdel = cos(del0);

  double salpsqr = s_alpha * s_alpha;
  double sbetasqr = s_beta * s_beta;  

  if(swz_3rrs_invkin2(alpha, beta, params, sol2)) 
    std::cout<<"Error: Check Inverse Kinematics";

  gamma2 = sol2[1];

  // Direction cosines of strut vector
  stVec2[0] = -cos(gamma2)/2.;
  stVec2[1] = cos(gamma2)*SQRT3/2.;
  stVec2[2] = sin(gamma2);

  double temp1 = 1 + c_alpha*c_beta;
  double temp = sqrt(temp1 * temp1);

  // Direction cosines of wedge vector
  vecWedge2[0] = cdel*s_beta + (c_beta*(c_alpha + c_beta)*sdel)/
    (2.*temp) - 
   (SQRT3*c_beta*s_alpha*s_beta*sdel)/
    (2.*temp);
  vecWedge2[1] = -(c_beta*cdel*s_alpha) + ((-((c_alpha*s_alpha*s_beta)/
           temp) + 
        ((c_alpha + c_beta)*s_alpha*s_beta)/temp)
       *sdel)/2. - (SQRT3*((c_alpha*(c_alpha + c_beta))/
         temp + 
        (salpsqr*sbetasqr)/temp)*
      sdel)/2.;
  vecWedge2[2] = c_alpha*c_beta*cdel + ((-((c_alpha*(c_alpha + c_beta)*s_beta)/
           temp) - 
        (salpsqr*s_beta)/temp)*sdel)/2.\
    - (SQRT3*(((c_alpha + c_beta)*s_alpha)/temp - 
        (c_alpha*s_alpha*sbetasqr)/temp)*
      sdel)/2.;   

  del2 = vecAngle(vecWedge2, stVec2);
  return 0;
}

/** swz_3rrs_invkin_delta3: Inverse Kinematics to find delta1, spherical joint angle of 3rd leg

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched
@param delta3 The variable to store the solution

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: 0, if solved properly

*/
int swz_3rrs_invkin_delta3(double alpha, double beta, double *params, double &del3)
{
  double del0;
  double sol3[2], gamma3;
  double vecWedge3[3], stVec3[3];

  del0 = params[5];

  double s_alpha = sin(alpha);
  double c_alpha = cos(alpha);
  double s_beta = sin(beta);
  double c_beta = cos(beta);

  double sdel = sin(del0);
  double cdel = cos(del0);

  double salpsqr = s_alpha * s_alpha;
  double sbetasqr = s_beta * s_beta;  

  if(swz_3rrs_invkin3(alpha, beta, params, sol3)) 
    std::cout<<"Error: Check Inverse Kinematics";

  gamma3 = sol3[1];

  // Direction cosines of strut vector
  stVec3[0] = -cos(gamma3)/2.;
  stVec3[1] = -cos(gamma3)*SQRT3/2.;
  stVec3[2] = sin(gamma3);

  double temp1 = 1 + c_alpha*c_beta;
  double temp = sqrt(temp1 * temp1);

  // Direction cosines of wedge vector
  vecWedge3[0] = cdel*s_beta + (c_beta*(c_alpha + c_beta)*sdel)/
    (2.*temp) + 
   (SQRT3*c_beta*s_alpha*s_beta*sdel)/
    (2.*temp);
  vecWedge3[1] = -(c_beta*cdel*s_alpha) + ((-((c_alpha*s_alpha*s_beta)/
           temp) + 
        ((c_alpha + c_beta)*s_alpha*s_beta)/temp)
       *sdel)/2. + (SQRT3*((c_alpha*(c_alpha + c_beta))/
         temp + 
        (salpsqr*sbetasqr)/temp)*
      sdel)/2.;
  vecWedge3[2] = c_alpha*c_beta*cdel + ((-((c_alpha*(c_alpha + c_beta)*s_beta)/
           temp) - 
        (salpsqr*s_beta)/temp)*sdel)/2.\
    + (SQRT3*(((c_alpha + c_beta)*s_alpha)/temp - 
        (c_alpha*s_alpha*sbetasqr)/temp)*
      sdel)/2.;   

  del3 = vecAngle(vecWedge3, stVec3);
  return 0;
}


/** S1 function for 1st leg: S11

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: value, value of the function

*/
double swz_3rrs_S11(double alpha, double beta, double *params)
{
	double theta1;

  if(swz_3rrs_invkin_theta1(alpha, beta, params, theta1))
    std::cout<<"Error: Check Inverse Kinematics";     
  return(theta1*theta1);

	// double rb,ra,z;
	// double value;
	// double theta1;
	
	// rb = params[0]; 
	// ra = params[3]; 
	// z = params[4];
	
 //  double c_alpha = cos(alpha), s_alpha = sin(alpha);
 //  double c_beta = cos(beta), s_beta = sin(beta);

 //  double calpsqr = c_alpha*c_alpha;
 //  double salpsqr = s_alpha * s_alpha;
 //  double cbetasqr =  c_beta * c_beta;
 //  double sbetasqr = s_beta * s_beta;

 //  if(swz_3rrs_invkin_theta1(alpha, beta, params, theta1))
 //    std::cout<<"Error: Check Inverse Kinematics";   
	
	// double s_theta1=sin(theta1), c_theta1=cos(theta1); 

	
	// value = 16*(1 + c_alpha*c_beta)*c_theta1*(z - ra*s_beta) + (-6*ra + 16*rb + \
	// 2*ra*calpsqr - 16*ra*c_alpha*c_beta + 16*rb*c_alpha*c_beta - \
	// 14*ra*cbetasqr + 2*ra*calpsqr*cbetasqr - \
	// 2*ra*salpsqr - 2*ra*cbetasqr*salpsqr + \
	// 14*ra*sbetasqr - 2*ra*calpsqr*sbetasqr + \
	// 2*ra*salpsqr*sbetasqr)*s_theta1;
	// //std::cout<<"S11 val: "<<value<<"\t alpha: "<<alpha<<"\t beta: "<<beta<<"\t z: "<<z<<std::endl;
	// return(value);
}


/** S1 function for 2nd leg: S12

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'

@returns: value, value of the function

*/
double swz_3rrs_S12(double alpha, double beta, double *params)
{
  double theta2;

  if(swz_3rrs_invkin_theta2(alpha, beta, params, theta2))
    std::cout<<"Error: Check Inverse Kinematics";     
  return(theta2*theta2);

 //  double rb,ra,z;
 //  double value;
 //  double theta2;
  
 //  rb = params[0]; 
 //  ra = params[3]; 
 //  z = params[4];
  
 //  double c_alpha = cos(alpha), s_alpha = sin(alpha);
 //  double c_beta = cos(beta), s_beta = sin(beta);

 //  double calpsqr = c_alpha*c_alpha;
 //  double salpsqr = s_alpha * s_alpha;
 //  double cbetasqr =  c_beta * c_beta;
 //  double sbetasqr = s_beta * s_beta;

 //  if(swz_3rrs_invkin_theta2(alpha, beta, params, theta2))
 //    std::cout<<"Error: Check Inverse Kinematics";   
  
 //  double s_theta2=sin(theta2), c_theta2=cos(theta2); 
	
	// value = -2*c_theta2*ra*s_beta + 3*ra*s_theta2 - 4*rb*s_theta2 - \
	// ra*salpsqr*s_theta2 - cbetasqr*ra*(1 + \
	// salpsqr)*s_theta2 + ra*sbetasqr*s_theta2 + \
	// ra*salpsqr*sbetasqr*s_theta2 + calpsqr*ra*(1 + \
	// cbetasqr - sbetasqr)*s_theta2 - \
	// 2*SQRT3*c_beta*ra*s_alpha*(c_theta2 + 2*s_beta*s_theta2) - \
	// 4*c_theta2*z - c_alpha*(4*c_beta*(-ra + rb)*s_theta2 + \
	// c_theta2*(ra*(2*c_beta*s_beta + SQRT3*s_alpha*(1 + cbetasqr - \
	// sbetasqr)) + 4*c_beta*z));
	// //std::cout<<"S11 val: "<<value<<"\t alpha: "<<alpha<<"\t beta: "<<beta<<"\t z: "<<z<<std::endl;
	// return(value);
}

/** S1 function for 3rd leg: S13

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'

@returns: value, value of the function

*/
double swz_3rrs_S13(double alpha,double beta,double *params)
{
  double theta3;

  if(swz_3rrs_invkin_theta3(alpha, beta, params, theta3))
    std::cout<<"Error: Check Inverse Kinematics";     
  return(theta3*theta3);

 //  double rb,ra,z;
 //  double value;
 //  double theta3;
  
 //  rb = params[0]; 
 //  ra = params[3]; 
 //  z = params[4];
  
 //  double c_alpha = cos(alpha), s_alpha = sin(alpha);
 //  double c_beta = cos(beta), s_beta = sin(beta);

 //  double calpsqr = c_alpha*c_alpha;
 //  double salpsqr = s_alpha * s_alpha;
 //  double cbetasqr =  c_beta * c_beta;
 //  double sbetasqr = s_beta * s_beta;

 //  if(swz_3rrs_invkin_theta3(alpha, beta, params, theta3))
 //    std::cout<<"Error: Check Inverse Kinematics"; 
	
	// double s_theta3=sin(theta3), c_theta3=cos(theta3);
	
	// value = -2*c_theta3*ra*s_beta + 3*ra*s_theta3 - 4*rb*s_theta3 - \
	// ra*salpsqr*s_theta3 - cbetasqr*ra*(1 + \
	// salpsqr)*s_theta3 + ra*sbetasqr*s_theta3 + \
	// ra*salpsqr*sbetasqr*s_theta3 + calpsqr*ra*(1 + \
	// cbetasqr - sbetasqr)*s_theta3 + \
	// 2*SQRT3*c_beta*ra*s_alpha*(c_theta3 + 2*s_beta*s_theta3) - \
	// 4*c_theta3*z + c_alpha*(4*c_beta*(ra - rb)*s_theta3 + \
	// c_theta3*(ra*(-2*c_beta*s_beta + SQRT3*s_alpha*(1 + cbetasqr - \
	// sbetasqr)) - 4*c_beta*z));
	// //std::cout<<"S11 val: "<<value<<"\t alpha: "<<alpha<<"\t beta: "<<beta<<"\t z: "<<z<<std::endl;
	// return(value);
}

/** S2 function: Gain singularity

@param alpha The angle representing pose of the end-effector plate
@param beta The angle representing pose of the end-effector plate
@param params The length data of the manipulator and the z-slice in which it is searched

Variables used:
Variable | Description
-------- |------------
b       | Location of the base of the crack 
l0       | Length of the coupler also distance between the bases of crank and rocker
l1       | Length of the crank
r        | Length of the strut 
rt       | Radius of the end-effector plate
lamda    | Ratio of the distance between the point on the coupler where the strut is located to the length of the coupler
del0     | The wedge angle of the axis of spherical joint with the the normal to end-effector plate
z        | Z-slice where the scanner searches in this iteration
qi       | Location of the end points of the traingular end-effector plate, i = 1, 2, 3
soli     | Solution for inverse kinematics containing theta and gamma for each leg, i = 1, 2, 3
s_i      | Sine of the variable 'i'
c_i      | Cosine of the variable 'i'
temp_i   | Temporary variable used for expressions evaluated repeatedly
det      | Value of the function, the determinant

@returns: det, value of the function

*/
double swz_3rrs_S2(double alpha, double beta, double *params)
{
	
	double rb,lcr,lst;
	double sol1[2], sol2[2], sol3[2];  
	
	rb = params[0]; 
	lcr = params[1]; 
	lst = params[2];
	
  if(swz_3rrs_invkin1(alpha, beta, params, sol1))
    std::cout<<"Error: Check SWZ";
  if(swz_3rrs_invkin2(alpha, beta, params, sol2))
    std::cout<<"Error: Check SWZ";  
  if(swz_3rrs_invkin3(alpha, beta, params, sol3))
    std::cout<<"Error: Check SWZ";

//gamma is the passive rotary joint angle
	double gamma1=sol1[1], gamma2=sol2[1], gamma3=sol3[1], theta1=sol1[0], theta2=sol2[0], theta3=sol3[0];
	
	double c_gamma1 = cos(gamma1), c_gamma2 = cos(gamma2), c_gamma3 = cos(gamma3);
	double s_gamma1 = sin(gamma1), s_gamma2 = sin(gamma2), s_gamma3 = sin(gamma3);
	
	double c_theta1 = cos(theta1), c_theta2 = cos(theta2), c_theta3 = cos(theta3);
	double s_theta1 = sin(theta1), s_theta2 = sin(theta2), s_theta3 = sin(theta3);
	
	double s_T2mG2 = sin(theta2-gamma2) , s_T1mG3 = sin(theta1-gamma3) , s_T3mG3 = sin(theta3-gamma3);
	double s_T1pG3 = sin(theta1+gamma3), s_G1pG3 = sin(gamma1+gamma3);
	
	
	double det = ((2*c_gamma2*(lcr*s_theta1 - lcr*s_theta2 + lst*s_gamma1) \
	+ (3*rb + lcr*c_theta1 + 2*lcr*c_theta2 + \
	lst*c_gamma1)*s_gamma2)*(2*c_gamma3*(lcr*s_theta2 - lcr*s_theta3 + \
	lst*s_gamma2) + (3*rb + lcr*c_theta2 + 2*lcr*c_theta3 + \
	lst*c_gamma2)*s_gamma3)*(-((3*rb + 2*lcr*c_theta1 + lcr*c_theta3 + \
	lst*c_gamma3)*s_gamma1) + 2*c_gamma1*(lcr*s_theta1 - lcr*s_theta3 - \
	lst*s_gamma3)) - (((3*rb + 2*lcr*c_theta1 + lcr*c_theta2 + \
	lst*c_gamma2)*s_gamma1 + c_gamma1*(-2*lcr*s_theta1 + 2*lcr*s_theta2 + \
	2*lst*s_gamma2))*(-2*lcr*s_T2mG2 + (3*rb + lcr*c_theta3 + \
	lst*c_gamma3)*s_gamma2 + 2*c_gamma2*(lcr*s_theta3 + \
	lst*s_gamma3))*(lcr*s_T1mG3 - 4*lcr*s_T3mG3 + lst*s_T1mG3 + \
	6*rb*s_gamma3 + 3*lcr*s_T1pG3 + 3*lst*s_G1pG3))/2.);
	
	
	
	return(det);
}

double swz_3rrs_S301(double alpha, double beta, double *params)
{
  double del1;
  if(swz_3rrs_invkin_delta1(alpha, beta, params, del1)) 
    std::cout<<"Check INVK";
  double val = DELTA_LIMIT_MAX - del1;

  if(val>0)
    return 1;
  else
    return 0;
}


double swz_3rrs_S302(double alpha, double beta, double *params)
{
  double del2;
  if(swz_3rrs_invkin_delta2(alpha, beta, params, del2)) 
    std::cout<<"Check INVK";
  double val = DELTA_LIMIT_MAX - del2;

  if(val>0)
    return 1;
  else
    return 0;
}


double swz_3rrs_S303(double alpha, double beta, double *params)
{
  double del3;
  if(swz_3rrs_invkin_delta3(alpha, beta, params, del3)) 
    std::cout<<"Check INVK";
  double val = DELTA_LIMIT_MAX - del3;

  if(val>0)
    return 1;
  else
    return 0;
}

