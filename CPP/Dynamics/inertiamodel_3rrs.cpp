#include "inertiamodel_3rrs.h"

InertiaModel::InertiaModel()
{
  g = 9.812;
  Eyng = 180*1e9; // in Pa
  density = 8050;//Density in kg/m3
  delSt = 0.1*1e-3;//Maximum allowable deflection of strut  
  delCr = 0.1*1e-3;//Maximum allowable deflection of crank
 
  h0Cr = (15+10)*1e-3;//Minimum Height of the crank
  h0St = (15+10)*1e-3;//Minimum Height of the strut

  minb = 30*1e-3; //Minimum value of breadth for the links

  tTp = 0.012;//Thickness of the top platform
  wTp = 0.080;//Width of the top platform

  sqrt3 = sqrt(3);
}

//Set the lengths from the inputs
void InertiaModel::setlength(double l[])
{
  b = l[0];
  lCr = l[1];
  lSt = l[2];
  rt = l[3];

  dt = rt*sqrt(3);  
}

void InertiaModel::setpayload(double mpl[3], double platform[3], double cylradius, double cylheight, double cyloffset, double cubside, double cubheight, double cuboffset)
{

  // This function accepts payloads and calculates total mass, center of mass and moment of inertia 
  double masspayload = mpl[0] + mpl[1] + mpl[2];
  
  P = (masspayload ) * g;

  double xPl1 = 0;
  double yPl1 = 0;
  double zPl1 = 0.1;

  double xPl2 = cyloffset;
  double yPl2 = 0;
  double zPl2 = cylheight/2 ;

  double xPl3 = cuboffset;
  double yPl3 = 0;
  double zPl3 = cubheight/2 ;
 
  // Inertia of Payload 1
  double IPl1x = mpl[0] * ((platform[1]*platform[1]+platform[2]*platform[2])/12 + (yPl1*yPl1 + zPl1*zPl1));
  double IPl1y = mpl[0] * ((platform[0]*platform[0]+platform[2]*platform[2])/12 + (xPl1*xPl1 + zPl1*zPl1));
  double IPl1z = mpl[0] * ((platform[1]*platform[1]+platform[0]*platform[0])/12 + (yPl1*yPl1 + xPl1*xPl1));

  // Inertia of Payload 2
  double IPl2x = mpl[1] * ((3*cylradius*cylradius + cylheight*cylheight)/12 + (yPl2*yPl2 + zPl2*zPl2));
  double IPl2y = mpl[1] * ((3*cylradius*cylradius + cylheight*cylheight)/12 + (xPl2*xPl2 + zPl2*zPl2));
  double IPl2z = mpl[1] * ((cylradius*cylradius/2) + (yPl2*yPl2 + xPl2*xPl2));  

  // Inertia of Payload 3
  double IPl3x = mpl[2] * ((cubside*cubside + cubheight*cubheight)/12 + (yPl3*yPl3 + zPl3*zPl3));
  double IPl3y = mpl[2] * ((cubside*cubside + cubheight*cubheight)/12 + (xPl3*xPl3 + zPl3*zPl3));
  double IPl3z = mpl[2] * ((cubside*cubside)/6 + (yPl3*yPl3 + xPl3*xPl3));

  inertia[3][0] = IPl1x + IPl2x + IPl3x;
  inertia[3][1] = IPl1y + IPl2y + IPl3y;  
  inertia[3][2] = IPl1z + IPl2z + IPl3z;

  mass[3] = masspayload; 

  xPl = (mpl[0]*xPl1 + mpl[1]*xPl2 + mpl[2]*xPl3)/masspayload;
  yPl = (mpl[0]*yPl1 + mpl[1]*yPl2 + mpl[2]*yPl3)/masspayload;
  zPl = (mpl[0]*zPl1 + mpl[1]*zPl2 + mpl[2]*zPl3)/masspayload;

//Compute the inertia matrices of crank, strut and moving triangle
  calc();

}

//Returns only the elbow out solution when the platform is above which would always be the case here
bool InertiaModel::rrdyad(double l1, double l2, double p[2][2], double soln[2][2])
{
  double del_x = p[1][0]-p[0][0]; 
  double del_y = p[1][1]-p[0][1];
  double d = del_y*del_y+del_x*del_x;
  double sqrtd = sqrt(d);
  if((sqrtd > l1 + l2) ||( sqrtd < std::abs(l1 - l2)))
    return 1;
  
  double alpha = atan2(del_y,del_x);
  double beta = acos((l1*l1-l2*l2+d)/(2.0*l1*sqrtd));
  soln[0][0] = alpha-beta;
  //	soln[1][0]=alpha+beta;
  //	soln[1]=atan2(del_y-l[0]*sin(soln[0]),del_x-l[0]*cos(soln[0]))-2*M_PI;	
  soln[0][1] = atan2(del_y-l1*sin(soln[0][0]),del_x-l1*cos(soln[0][0]));	
  //	soln[1][1]=atan2(-del_y+l1*sin(soln[1][0]),-del_x+l1*cos(soln[1][0]));
  return 0;
}

bool InertiaModel::invk(double task[3], double theta[3], double gamma[3])
{
 
  double c1 = task[0];
  double c2 = task[1];
  double z = task[2];

  double temp1 = c1*c1 - c2*c2;
  double delc = 1+c1*c1+c2*c2;
  double temp2 = rt/delc;

  double p1[3] = {(1+2*temp1)*temp2,0,-2*c2*temp2+z};
  double p2[3] = {(1-temp1-2*sqrt3*c1*c2)*temp2,0,((sqrt3*c1+c2)*temp2+z)};
  double p3[3] = {(1-temp1+2*sqrt3*c1*c2)*temp2,0,((-sqrt3*c1+c2)*temp2+z)};

  double sol[2][2]={{0,0},{0,0}};
  double rrpnt1[2][2]={{b,0},{p1[0],p1[2]}};
  if(rrdyad(lCr, lSt, rrpnt1, sol))
    return 1;
  theta[0] = sol[0][0];
  gamma[0] = sol[0][1];

  double rrpnt2[2][2]={{b,0},{p2[0],p2[2]}};
  if(rrdyad(lCr, lSt, rrpnt2, sol))
    return 1;
  theta[1] = sol[0][0];
  gamma[1] = sol[0][1];

  double rrpnt3[2][2]={{b,0},{p3[0],p3[2]}};
  if(rrdyad(lCr, lSt, rrpnt3, sol))
    return 1;
  theta[2] = sol[0][0];
  gamma[2] = sol[0][1];

  return 0;
  
}

//Firs order and second order inverse kinematics:: the joint rates and joint accelerations would be returned
bool InertiaModel::invk1(double phi[3], double omega[3], double alpha[3], double zderv[3], VectorXd &q9, VectorXd &dq9, VectorXd &ddq9)
{
  double c1, c2, z;
  double tanphiby2 = tan(phi[2]/2.);
  c1 = tanphiby2*phi[0];
  c2 = tanphiby2*phi[1];
  z = zderv[0];

  double theta[3];
  double gamma[3];
  double task[3] =
    {
      c1, c2, z
    }
  ;
  
  if(  invk(task, theta, gamma))
    {
      std::cout << "Error: Check SWZ"<< std::endl;
      return 1;
    }
  
  double dc1, dc2, dz;
  double consdtask = 1./(cos(phi[2]/2.)*cos(phi[2]/2.)*2.);
  dc1 = consdtask*omega[0]*omega[2];
  dc2 = consdtask*omega[1]*omega[2];
  dz = zderv[1];
  
  double ddc1, ddc2, ddz;
  double consddtask = consdtask*(tanphiby2*omega[2]*omega[2]+alpha[2]);
  ddc1 = consddtask*alpha[0];
  ddc2 = consddtask*alpha[1];
  ddz = zderv[2];
  
//  double lCrsqr, lStsqr, rtsqr, c1sqr, c2sqr, delc, delcsqr;

//  lCrsqr = lCr*lCr;
//  lStsqr = lSt*lSt;
//  rtsqr = rt*rt;
  double c1sqr = c1*c1;
  double c2sqr = c2*c2;

  double delc = 1+c1sqr + c2sqr;
//  delcsqr = delc*delc;

  double tht1 = theta[0];
  double tht2 = theta[1];
  double tht3 = theta[2];

  double gam1 = gamma[0];
  double gam2 = gamma[1];
  double gam3 = gamma[2];
  
  double stht1, stht2, stht3, ctht1, ctht2, ctht3, sgam1, sgam2, sgam3, cgam1, cgam2, cgam3;
  
  stht1 = sin(tht1);
  stht2 = sin(tht2);
  stht3 = sin(tht3);

  ctht1 = cos(tht1);
  ctht2 = cos(tht2);
  ctht3 = cos(tht3);

  sgam1 = sin(gam1);
  sgam2 = sin(gam2);
  sgam3 = sin(gam3);

  cgam1 = cos(gam1);
  cgam2 = cos(gam2);
  cgam3 = cos(gam3);

  
  MatrixXd JetaBeta(6,6);

  JetaBeta << delc*lCr*stht1, 0, 0, delc*lSt*sgam1, 0, 0, 
              (-ctht1)*delc*lCr, 0, 0, (-cgam1)*delc*lSt, 0, 0, 
              0, delc*lCr*stht2, 0, 0, delc*lSt*sgam2, 0, 
              0, (-ctht2)*delc*lCr, 0, 0, (-cgam2)*delc*lSt, 0, 
              0, 0, delc*lCr*stht3, 0, 0, delc*lSt*sgam3, 
              0, 0, (-ctht3)*delc*lCr, 0, 0, (-cgam3)*delc*lSt;
     
  MatrixXd JetaX(6,3);
   
  JetaX <<  -2*c1*ctht1*lCr - 2*c1*cgam1*lSt + 4*c1*rt -  2*c1*b, -2*c2*ctht1*lCr - 2*c2*cgam1*lSt - 4*c2*rt - 2*c2*b, 0, 
            -2*c1*lSt*sgam1 - 2*c1*lCr*stht1 + 2*c1*z, -2*rt - 2*c2*lSt*sgam1 - 2*c2*lCr*stht1 + 2*c2*z, delc, 
            -2*c1*ctht2*lCr - 2*c1*cgam2*lSt - 2*c1*rt - 2*c1*b - 2*c2*rt*sqrt3, -2*c2*ctht2*lCr - 2*c2*cgam2*lSt + 2*c2*rt - 2*c2*b - 2*c1*rt*sqrt3, 0, 
            -2*c1*lSt*sgam2 + rt*sqrt3 - 2*c1*lCr*stht2 + 2*c1*z, rt - 2*c2*lSt*sgam2 - 2*c2*lCr*stht2 + 2*c2*z, delc, 
            -2*c1*ctht3*lCr - 2*c1*cgam3*lSt - 2*c1*rt - 2*c1*b + 2*c2*rt*sqrt3, -2*c2*ctht3*lCr - 2*c2*cgam3*lSt + 2*c2*rt - 2*c2*b + 2*c1*rt*sqrt3, 0, 
            -2*c1*lSt*sgam3 - rt*sqrt3 - 2*c1*lCr*stht3 + 2*c1*z, rt - 2*c2*lSt*sgam3 - 2*c2*lCr*stht3 + 2*c2*z, delc;


   Vector3d dX;
   dX <<  dc1, dc2, dz;
   
   VectorXd dBeta(6);
   dBeta = -JetaBeta.inverse()*(JetaX*dX);

   //   VectorXd q9(9);
   q9 << theta[0], theta[1], theta[2], gamma[0], gamma[1], gamma[2],  c1, c2, z;

   // VectorXd dq9(9);
   dq9 << dBeta, dc1, dc2, dz;
   

   double dtht1, dtht2, dtht3, dgam1, dgam2, dgam3;
   dtht1 = dq9(0);
   dtht2 = dq9(1);
   dtht3 = dq9(2);

   dgam1 = dq9(3);
   dgam2 = dq9(4);
   dgam3 = dq9(5);
   
   MatrixXd dJetaBeta(6,6);

   dJetaBeta << ctht1*delc*dtht1*lCr + (2*c1*dc1 + 2*c2*dc2)*lCr*stht1, 0, 0, cgam1*delc*dgam1*lSt + (2*c1*dc1 + 2*c2*dc2)*lSt*sgam1, 0, 0, 
                (-ctht1)*(2*c1*dc1 + 2*c2*dc2)*lCr + delc*dtht1*lCr*stht1, 0, 0, (-cgam1)*(2*c1*dc1 + 2*c2*dc2)*lSt + delc*dgam1*lSt*sgam1, 0, 0, 
                0, ctht2*delc*dtht2*lCr + (2*c1*dc1 + 2*c2*dc2)*lCr*stht2, 0, 0, cgam2*delc*dgam2*lSt + (2*c1*dc1 + 2*c2*dc2)*lSt*sgam2, 0, 
                0, (-ctht2)*(2*c1*dc1 + 2*c2*dc2)*lCr + delc*dtht2*lCr*stht2, 0, 0, (-cgam2)*(2*c1*dc1 + 2*c2*dc2)*lSt + delc*dgam2*lSt*sgam2, 0, 
                0, 0, ctht3*delc*dtht3*lCr + (2*c1*dc1 + 2*c2*dc2)*lCr*stht3, 0, 0, cgam3*delc*dgam3*lSt + (2*c1*dc1 + 2*c2*dc2)*lSt*sgam3, 
                0, 0, (-ctht3)*(2*c1*dc1 + 2*c2*dc2)*lCr + delc*dtht3*lCr*stht3, 0, 0, (-cgam3)*(2*c1*dc1 + 2*c2*dc2)*lSt + delc*dgam3*lSt*sgam3;

   Vector3d ddX;
   ddX <<  ddc1, ddc2, ddz;

   MatrixXd dJetaX(6,3);

    dJetaX <<   -2*ctht1*dc1*lCr - 2*cgam1*dc1*lSt + 4*dc1*rt - 2*dc1*b + 2*c1*dgam1*lSt*sgam1 + 2*c1*dtht1*lCr*stht1, 
   -2*ctht1*dc2*lCr - 2*cgam1*dc2*lSt - 4*dc2*rt - 2*dc2*b + 2*c2*dgam1*lSt*sgam1 + 2*c2*dtht1*lCr*stht1, 0 , 
   2*c1*dz - 2*c1*ctht1*dtht1*lCr - 2*c1*cgam1*dgam1*lSt - 2*dc1*lSt*sgam1 - 2*dc1*lCr*stht1 + 2*dc1*z, 
   2*c2*dz - 2*c2*ctht1*dtht1*lCr - 2*c2*cgam1*dgam1*lSt - 2*dc2*lSt*sgam1 - 2*dc2*lCr*stht1 + 2*dc2*z, 2*c1*dc1 + 2*c2*dc2 , 
   -2*ctht2*dc1*lCr - 2*cgam2*dc1*lSt - 2*dc1*rt - 2*sqrt3*dc2*rt - 2*dc1*b + 2*c1*dgam2*lSt*sgam2 + 2*c1*dtht2*lCr*stht2, 
   -2*ctht2*dc2*lCr - 2*cgam2*dc2*lSt - 2*sqrt3*dc1*rt + 2*dc2*rt - 2*dc2*b + 2*c2*dgam2*lSt*sgam2 + 2*c2*dtht2*lCr*stht2, 0 , 
   2*c1*dz - 2*c1*ctht2*dtht2*lCr - 2*c1*cgam2*dgam2*lSt - 2*dc1*lSt*sgam2 - 2*dc1*lCr*stht2 + 2*dc1*z, 
   2*c2*dz - 2*c2*ctht2*dtht2*lCr - 2*c2*cgam2*dgam2*lSt - 2*dc2*lSt*sgam2 - 2*dc2*lCr*stht2 + 2*dc2*z, 2*c1*dc1 + 2*c2*dc2 , 
   -2*ctht3*dc1*lCr - 2*cgam3*dc1*lSt - 2*dc1*rt + 2*sqrt3*dc2*rt - 2*dc1*b + 2*c1*dgam3*lSt*sgam3 + 2*c1*dtht3*lCr*stht3, 
   -2*ctht3*dc2*lCr - 2*cgam3*dc2*lSt + 2*sqrt3*dc1*rt + 2*dc2*rt - 2*dc2*b + 2*c2*dgam3*lSt*sgam3 + 2*c2*dtht3*lCr*stht3, 0 , 
   2*c1*dz - 2*c1*ctht3*dtht3*lCr - 2*c1*cgam3*dgam3*lSt - 2*dc1*lSt*sgam3 - 2*dc1*lCr*stht3 + 2*dc1*z, 
   2*c2*dz - 2*c2*ctht3*dtht3*lCr - 2*c2*cgam3*dgam3*lSt - 2*dc2*lSt*sgam3 - 2*dc2*lCr*stht3 + 2*dc2*z, 2*c1*dc1 + 2*c2*dc2 ;

   VectorXd ddBeta(6);
   ddBeta = -JetaBeta.inverse()*(JetaX*ddX + dJetaX*dX + dJetaBeta*dBeta);

   // VectorXd ddq9(9);
   ddq9 <<  ddBeta, ddc1, ddc2, ddz;

   return 0;

}

void InertiaModel::calc()
{
  
  double hCr, hSt;//Height of the corresponding links
  double bCr, bSt;//Width of the corresponding links

//Set the breadth of crank and strut to be minimum values  
  bCr = minb;
  bSt = minb;
 
  //Find the height for stipulated deflection
  hCr = lCr*std::pow(4.0*P/(bCr*Eyng*delCr),1/3.);

  //If the obtained height is less than min allowable height, then use min height for further calculations 
  if(hCr<h0Cr)
    hCr = h0Cr;

//Inertia matrix computation for the crank
  mass[0] = bCr*hCr*lCr*density; // Crank mass

  inertia[0][0] = mass[0]*(hCr*hCr + bCr*bCr)/12.;
  inertia[0][1] = mass[0]*(lCr*lCr + hCr*hCr)/12.;
  inertia[0][2] = mass[0]*(lCr*lCr + bCr*bCr)/12.;

//Inertia matrix computation for the strut
  hSt= lSt*std::pow(4.0*P/(bSt*Eyng*delSt),1/3.);
  if(hSt<h0St)
   hSt = h0St;
  
  mass[1] = bSt*hSt*lSt*density; // Strut mass
  inertia[1][0] = mass[1]*(bSt*bSt + hSt*hSt)/12.;
  inertia[1][1] = mass[1]*(lSt*lSt + hSt*hSt)/12.;
  inertia[1][2] = mass[1]*(lSt*lSt + bSt*bSt)/12.;

//Inertia matrix computation for the top platform
  mass[2] = 3*sqrt3*rt*tTp*wTp*density; //Top plate mass

  inertia[2][0] = mass[2]*(3*((rt- wTp/2)*(rt - wTp/2)+(rt + wTp/2)*(rt + wTp/2))+tTp*tTp)/12;
  inertia[2][1] = inertia[2][0];
  inertia[2][2] = mass[2]*((rt- wTp/2)*(rt - wTp/2)+(rt + wTp/2)*(rt + wTp/2))/2;
}
void InertiaModel::output()
{
  std::cout<<"\nb: "<<b<<"\nlCr: "<<lCr<<"\nlSt: "<<lSt<<"\nrt: "<<rt<<std::endl;
  std::cout<<"\nmass 1: "<<mass[0]<<"\nmass 2: "<<mass[1]<<"\nmass 3: "<<mass[2]<<"\nmass 4: "<<mass[3] ;
  std::cout<<"\nICrx: "<<inertia[0][0]<<"\nICry: "<<inertia[0][1]<<"\nICrz: "<<inertia[0][2];
  std::cout<<"\nIStx: "<<inertia[1][0]<<"\nISty: "<<inertia[1][1]<<"\nIStz: "<<inertia[1][2];
  std::cout<<"\nITpx: "<<inertia[2][0]<<"\nITpy: "<<inertia[2][1]<<"\nITpz: "<<inertia[2][2];
  std::cout<<"\nIPlx: "<<inertia[3][0]<<"\nIPly: "<<inertia[3][1]<<"\nIPlz: "<<inertia[3][2];
  std::cout<<"\nxPl: "<<xPl<<"\nyPl: "<<yPl<<"\nyPl: "<<zPl<<std::endl;
}
