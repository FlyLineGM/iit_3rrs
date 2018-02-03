#include "inertiamodel_3rrs.h"

void InertiaModel::Mq(VectorXd &q9, MatrixXd &MqMat)
{

  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  

  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double c1 = q9(6);
  double c2 = q9(7);


  double ICry = inertia[0][1];
  double ISty = inertia[1][1];

  double ITpx = inertia[2][0];
  double ITpy = inertia[2][1];
  double ITpz = inertia[2][2];

  double IPlx = inertia[3][0];
  double IPly = inertia[3][1];
  double IPlz = inertia[3][2];

  double mCr = mass[0];
  double mSt = mass[1];
  double mTp = mass[2];
  double mPl = mass[3];


  double lCrsqr = lCr*lCr;
  double lStsqr = lSt*lSt;
  double rtsqr = rt*rt;
  double c1sqr = c1*c1;
  double c2sqr = c2*c2;
  double c1pow4 = c1sqr*c1sqr;
  double c2pow4 = c2sqr*c2sqr;

  double delc = 1+c1sqr + c2sqr;
  double delcsqr = delc*delc;

  double cthphi1 = cos(tht1 - phi1);
  double cthphi2 = cos(tht2 - phi2);
  double cthphi3 = cos(tht3 - phi3);

  double xPlsqr = xPl*xPl;
  double yPlsqr = yPl*yPl;
  double zPlsqr = zPl*zPl;

  MqMat <<   
          ICry + (lCrsqr*(mCr + 4*mSt))/4.,0,0,(lCr*lSt*mSt*cthphi1)/2.,0,0,0,0, 0,
          0,ICry + (lCrsqr*(mCr + 4*mSt))/4.,0,0,(lCr*lSt*mSt*cthphi2)/2.,0,0,0, 0,
          0,0,ICry + (lCrsqr*(mCr + 4*mSt))/4.,0,0,(lCr*lSt*mSt*cthphi3)/2.,0,0,0,
          (lCr*lSt*mSt*cthphi1)/2.,0,0,ISty + (lStsqr*mSt)/4.,0,0,0,0,0,
          0,(lCr*lSt*mSt*cthphi2)/2.,0,0,ISty + (lStsqr*mSt)/4.,0,0,0,0,
          0,0,(lCr*lSt*mSt*cthphi3)/2.,0,0,ISty + (lStsqr*mSt)/4.,0,0,0,


          0,0,0,0,0,0,(4*(delcsqr*IPlx + ITpx + 2*c1sqr*ITpx + c1pow4*ITpx + \
          c1sqr*mPl*rtsqr + c1sqr*mTp*rtsqr + mPl*yPlsqr + 2*c1sqr*mPl*yPlsqr + \
          c1pow4*mPl*yPlsqr + pow(c2,6)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - \
          2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + mPl*zPlsqr + \
          2*c1sqr*mPl*zPlsqr + c1pow4*mPl*zPlsqr + 2*pow(c2,5)*mPl*(4*c1*rt*yPl \
          + (rt - xPl)*zPl) + 2*c2*mPl*(3*c1*rt*yPl - 3*pow(c1,3)*rt*yPl + (rt \
          - xPl)*zPl + c1pow4*(rt - xPl)*zPl - 2*c1sqr*(2*rt + xPl)*zPl) - \
          2*pow(c2,3)*mPl*(-7*c1*rt*yPl + 4*pow(c1,3)*rt*yPl + 2*(-rt + \
          xPl)*zPl + 2*c1sqr*(3*rt + xPl)*zPl) + c2pow4*(2*(1 + c1sqr)*IPlz + \
          ITpx + 2*ITpz + 2*c1sqr*ITpz + 2*mPl*rtsqr + 2*c1sqr*mPl*rtsqr + \
          2*mTp*rtsqr + 2*c1sqr*mTp*rtsqr - 4*mPl*rt*xPl + 12*c1sqr*mPl*rt*xPl \
          + 2*mPl*xPlsqr + 2*c1sqr*mPl*xPlsqr + 3*mPl*yPlsqr + \
          2*c1sqr*mPl*yPlsqr + mPl*zPlsqr) + c2sqr*(pow(1 + c1sqr,2)*IPlz + \
          2*(1 + c1sqr)*ITpx + ITpz + 2*c1sqr*ITpz + c1pow4*ITpz + mPl*rtsqr + \
          2*c1sqr*mPl*rtsqr + c1pow4*mPl*rtsqr + mTp*rtsqr + 2*c1sqr*mTp*rtsqr \
          + c1pow4*mTp*rtsqr - 2*mPl*rt*xPl + 8*c1sqr*mPl*rt*xPl - \
          2*c1pow4*mPl*rt*xPl + mPl*xPlsqr + 2*c1sqr*mPl*xPlsqr + \
          c1pow4*mPl*xPlsqr + 3*mPl*yPlsqr + 4*c1sqr*mPl*yPlsqr + \
          c1pow4*mPl*yPlsqr + 2*mPl*zPlsqr + \
          2*c1sqr*mPl*zPlsqr)))/pow(delc,4),(-4*(c1pow4*mPl*yPl*(-((3 + \
          8*c2sqr)*rt) + xPl + c2*zPl) + (1 + c2sqr)*mPl*yPl*(xPl + c2sqr*(rt + \
          xPl) + c2*zPl + pow(c2,3)*zPl) + pow(c1,5)*(c2*(IPlz + ITpz + \
          mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
          mPl*(rt - xPl)*zPl) + pow(c1,3)*(2*c2*(IPlz + ITpz + mPl*rtsqr + \
          mTp*rtsqr + mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 2*pow(c2,3)*(IPlz \
          + ITpz + mPl*rtsqr + mTp*rtsqr + 6*mPl*rt*xPl + mPl*xPlsqr + \
          mPl*yPlsqr) - 2*c2sqr*mPl*(5*rt + xPl)*zPl - mPl*(rt + 2*xPl)*zPl) + \
          c1*(pow(c2,5)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + \
          mPl*xPlsqr + mPl*yPlsqr) + 2*pow(c2,3)*(IPlz + ITpz + mPl*rtsqr + \
          mTp*rtsqr + 3*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + c2*(IPlz + ITpz \
          + mPl*(4*rt*xPl + xPlsqr + yPlsqr)) + c2pow4*mPl*(5*rt - xPl)*zPl - \
          mPl*(2*rt + xPl)*zPl - c2sqr*mPl*(rt + 2*xPl)*zPl) + \
          c1sqr*mPl*yPl*((-3 + 2*c2sqr + 8*c2pow4)*rt + 2*(1 + c2sqr)*(xPl + \
          c2*zPl))))/pow(delc,4),(-2*mPl*(c1sqr*yPl - (1 + c2sqr)*yPl + \
          c1*(-2*c2*xPl + 2*zPl)))/delcsqr,


          0,0,0,0,0,0,(-4*(c1pow4*mPl*yPl*(-((3 + 8*c2sqr)*rt) + xPl + c2*zPl) \
          + (1 + c2sqr)*mPl*yPl*(xPl + c2sqr*(rt + xPl) + c2*zPl + \
          pow(c2,3)*zPl) + pow(c1,5)*(c2*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - \
          2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + mPl*(rt - xPl)*zPl) + \
          pow(c1,3)*(2*c2*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + mPl*rt*xPl + \
          mPl*xPlsqr + mPl*yPlsqr) + 2*pow(c2,3)*(IPlz + ITpz + mPl*rtsqr + \
          mTp*rtsqr + 6*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) - \
          2*c2sqr*mPl*(5*rt + xPl)*zPl - mPl*(rt + 2*xPl)*zPl) + \
          c1*(pow(c2,5)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + \
          mPl*xPlsqr + mPl*yPlsqr) + 2*pow(c2,3)*(IPlz + ITpz + mPl*rtsqr + \
          mTp*rtsqr + 3*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + c2*(IPlz + ITpz \
          + mPl*(4*rt*xPl + xPlsqr + yPlsqr)) + c2pow4*mPl*(5*rt - xPl)*zPl - \
          mPl*(2*rt + xPl)*zPl - c2sqr*mPl*(rt + 2*xPl)*zPl) + \
          c1sqr*mPl*yPl*((-3 + 2*c2sqr + 8*c2pow4)*rt + 2*(1 + c2sqr)*(xPl + \
          c2*zPl))))/pow(delc,4),(4*(delcsqr*IPly + ITpy + 2*c2sqr*ITpy + \
          c2pow4*ITpy + c2sqr*mPl*rtsqr + c2sqr*mTp*rtsqr + 4*c2sqr*mPl*rt*xPl \
          + mPl*xPlsqr + 2*c2sqr*mPl*xPlsqr + c2pow4*mPl*xPlsqr + \
          pow(c1,6)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + \
          mPl*xPlsqr + mPl*yPlsqr) - 2*c2*mPl*rt*zPl + 2*pow(c2,3)*mPl*rt*zPl + \
          mPl*zPlsqr + 2*c2sqr*mPl*zPlsqr + c2pow4*mPl*zPlsqr + \
          2*pow(c1,5)*mPl*yPl*(-4*c2*rt + zPl) + 2*pow(c1,3)*mPl*yPl*(-5*c2*rt \
          + 4*pow(c2,3)*rt + 2*zPl + 2*c2sqr*zPl) + 2*c1*mPl*yPl*(-(c2*rt) + \
          pow(c2,3)*rt + zPl + 2*c2sqr*zPl + c2pow4*zPl) + c1pow4*(2*(1 + \
          c2sqr)*IPlz + ITpy + 2*ITpz + 2*c2sqr*ITpz + 2*mPl*rtsqr + \
          2*c2sqr*mPl*rtsqr + 2*mTp*rtsqr + 2*c2sqr*mTp*rtsqr - 4*mPl*rt*xPl + \
          12*c2sqr*mPl*rt*xPl + 3*mPl*xPlsqr + 2*c2sqr*mPl*xPlsqr + \
          2*mPl*yPlsqr + 2*c2sqr*mPl*yPlsqr - 8*c2*mPl*rt*zPl + mPl*zPlsqr) + \
          c1sqr*(pow(1 + c2sqr,2)*IPlz + 2*(1 + c2sqr)*ITpy + ITpz + \
          2*c2sqr*ITpz + c2pow4*ITpz + mPl*rtsqr + 2*c2sqr*mPl*rtsqr + \
          c2pow4*mPl*rtsqr + mTp*rtsqr + 2*c2sqr*mTp*rtsqr + c2pow4*mTp*rtsqr - \
          2*mPl*rt*xPl + 16*c2sqr*mPl*rt*xPl - 2*c2pow4*mPl*rt*xPl + \
          3*mPl*xPlsqr + 4*c2sqr*mPl*xPlsqr + c2pow4*mPl*xPlsqr + mPl*yPlsqr + \
          2*c2sqr*mPl*yPlsqr + c2pow4*mPl*yPlsqr - 10*c2*mPl*rt*zPl + \
          8*pow(c2,3)*mPl*rt*zPl + 2*mPl*zPlsqr + \
          2*c2sqr*mPl*zPlsqr)))/pow(delc,4),(-2*mPl*((1 + c1sqr - c2sqr)*xPl + \
          2*c2*(c1*yPl + zPl)))/delcsqr,


          0,0,0,0,0,0,(-2*mPl*(c1sqr*yPl - (1 + c2sqr)*yPl + c1*(-2*c2*xPl + \
          2*zPl)))/delcsqr,(-2*mPl*((1 + c1sqr - c2sqr)*xPl + 2*c2*(c1*yPl + \
          zPl)))/delcsqr,mPl + mTp;

}

void InertiaModel::Cq(VectorXd &q9, VectorXd &dq9, MatrixXd &CqMat)
{

//State variables
  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  

  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double c1 = q9(6);
  double c2 = q9(7);


// Derivatives of state variables
  double dtht1 = dq9(0);
  double dtht2 = dq9(1);
  double dtht3 = dq9(2);
  
  double dphi1 = dq9(3);
  double dphi2 = dq9(4);
  double dphi3 = dq9(5);

  double dc1 = dq9(6);
  double dc2 = dq9(7);

//Mass parameters
 
  double ITpx = inertia[2][0];
  double ITpy = inertia[2][1];
  double ITpz = inertia[2][2];

  double IPlx = inertia[3][0];
  double IPly = inertia[3][1];
  double IPlz = inertia[3][2];

  double mSt = mass[1];
  double mTp = mass[2];
  double mPl = mass[3];

//Simplification substitutions
  double rtsqr = rt*rt;
  double c1sqr = c1*c1;
  double c2sqr = c2*c2;
  double c1cub = c1*c1*c1;
  double c2cub = c2*c2*c2;
  double c1pow4 = c1sqr*c1sqr;
  double c2pow4 = c2sqr*c2sqr;

  double delc = 1+c1sqr + c2sqr;
 
  double sthphi1 = sin(tht1 - phi1);
  double sthphi2 = sin(tht2 - phi2);
  double sthphi3 = sin(tht3 - phi3);

  double xPlsqr = xPl*xPl;
  double yPlsqr = yPl*yPl;
  double zPlsqr = zPl*zPl;

  CqMat << 

        0,0,0,(dphi1*lCr*lSt*mSt*sthphi1)/2.,0,0,0,0,0,

        0,0,0,0,(dphi2*lCr*lSt*mSt*sthphi2)/2.,0,0,0,0,

        0,0,0,0,0,(dphi3*lCr*lSt*mSt*sthphi3)/2.,0,0,0,

        -(dtht1*lCr*lSt*mSt*sthphi1)/2.,0,0,0,0,0,0,0,0,

        0,-(dtht2*lCr*lSt*mSt*sthphi2)/2.,0,0,0,0,0,0,0,

        0,0,-(dtht3*lCr*lSt*mSt*sthphi3)/2.,0,0,0,0,0,0,

        0,0,0,0,0,0,(-4*(-(pow(c1,6)*dc2*(c2*(IPlz + ITpz + mPl*rtsqr + \
        mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + mPl*(rt - \
        xPl)*zPl)) + pow(c1,5)*(3*(1 + 4*c2sqr)*dc2*mPl*rt*yPl + 2*dc1*(IPlx \
        + ITpx + mPl*yPlsqr + c2sqr*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - \
        2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 2*c2*mPl*(rt - xPl)*zPl + \
        mPl*zPlsqr)) + pow(1 + c2sqr,2)*(c2cub*(-4*dc1*mPl*rt*yPl + \
        dc2*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr \
        + mPl*yPlsqr)) + 3*c2sqr*dc2*mPl*(rt - xPl)*zPl + dc2*mPl*(-rt + \
        xPl)*zPl - c2*(3*dc1*mPl*rt*yPl + dc2*(-2*IPlx + IPlz - 2*ITpx + ITpz \
        + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr - mPl*yPlsqr - \
        2*mPl*zPlsqr))) - c1pow4*(c2cub*(20*dc1*mPl*rt*yPl + dc2*(IPlz + \
        ITpz + mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + mPl*xPlsqr + \
        mPl*yPlsqr)) + c2sqr*dc2*mPl*(-25*rt + xPl)*zPl - 3*dc2*mPl*(rt + \
        xPl)*zPl + c2*(15*dc1*mPl*rt*yPl + dc2*(-2*IPlx + 3*IPlz - 2*ITpx + \
        3*ITpz + 3*mPl*rtsqr + 3*mTp*rtsqr + 6*mPl*rt*xPl + 3*mPl*xPlsqr + \
        mPl*yPlsqr - 2*mPl*zPlsqr))) + c1*(1 + c2sqr)*(3*(-1 + c2sqr + \
        4*c2pow4)*dc2*mPl*rt*yPl + dc1*(2*IPlx + 2*ITpx - mPl*rtsqr - \
        mTp*rtsqr + 2*mPl*yPlsqr + 2*c2pow4*(IPlz + ITpz + mPl*rtsqr + \
        mTp*rtsqr - 10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 4*c2*mPl*(4*rt \
        - xPl)*zPl + 4*c2cub*mPl*(5*rt - xPl)*zPl + 2*mPl*zPlsqr + \
        2*c2sqr*(IPlx + IPlz + ITpx + ITpz + mPl*rtsqr + mTp*rtsqr - \
        8*mPl*rt*xPl + mPl*xPlsqr + 2*mPl*yPlsqr + mPl*zPlsqr))) + \
        c1cub*(-10*c2sqr*(3 + 4*c2sqr)*dc2*mPl*rt*yPl + dc1*(4*IPlx + \
        4*ITpx + 3*mPl*rtsqr + 3*mTp*rtsqr + 4*mPl*yPlsqr + 4*c2pow4*(IPlz + \
        ITpz + mPl*rtsqr + mTp*rtsqr + 10*mPl*rt*xPl + mPl*xPlsqr + \
        mPl*yPlsqr) - 8*c2cub*mPl*(5*rt + xPl)*zPl - 4*c2*mPl*(7*rt + \
        2*xPl)*zPl + 4*mPl*zPlsqr + 4*c2sqr*(IPlx + IPlz + ITpx + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 7*mPl*rt*xPl + mPl*xPlsqr + 2*mPl*yPlsqr + \
        mPl*zPlsqr))) + c1sqr*(pow(c2,5)*(40*dc1*mPl*rt*yPl + dc2*(IPlz + \
        ITpz + mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + mPl*xPlsqr + \
        mPl*yPlsqr)) + 3*dc2*mPl*(rt + xPl)*zPl - 5*c2pow4*dc2*mPl*(7*rt + \
        xPl)*zPl - 2*c2sqr*dc2*mPl*(8*rt + xPl)*zPl - \
        2*c2cub*(-35*dc1*mPl*rt*yPl + dc2*(-2*IPlx + IPlz - 2*ITpx + ITpz \
        + mPl*rtsqr + mTp*rtsqr - 4*mPl*rt*xPl + mPl*xPlsqr - mPl*yPlsqr - \
        2*mPl*zPlsqr)) + c2*(30*dc1*mPl*rt*yPl + dc2*(4*IPlx - 3*IPlz + \
        4*ITpx - 3*ITpz + mPl*rtsqr + mTp*rtsqr - 6*mPl*rt*xPl - 3*mPl*xPlsqr \
        + mPl*yPlsqr + 4*mPl*zPlsqr)))))/pow(delc,5),(4*(dc1*(-3*pow(c1,5)*(1 \
        + 4*c2sqr)*mPl*rt*yPl + 10*c1cub*c2sqr*(3 + 4*c2sqr)*mPl*rt*yPl - \
        3*c1*(-1 + 5*c2pow4 + 4*pow(c2,6))*mPl*rt*yPl + pow(c1,6)*(c2*(IPlz + \
        ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + \
        mPl*yPlsqr) + mPl*(rt - xPl)*zPl) + c1pow4*(c2cub*(IPlz + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        c2sqr*mPl*(-25*rt + xPl)*zPl - 3*mPl*(rt + xPl)*zPl + c2*(-2*IPlx + \
        3*IPlz - 2*ITpx + 3*ITpz + 3*mPl*rtsqr + 3*mTp*rtsqr + 6*mPl*rt*xPl + \
        3*mPl*xPlsqr + mPl*yPlsqr - 2*mPl*zPlsqr)) - pow(1 + \
        c2sqr,2)*(c2cub*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - \
        2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 3*c2sqr*mPl*(rt - xPl)*zPl \
        + mPl*(-rt + xPl)*zPl + c2*(2*IPlx - IPlz + 2*ITpx - ITpz - mPl*rtsqr \
        - mTp*rtsqr + 2*mPl*rt*xPl - mPl*xPlsqr + mPl*yPlsqr + 2*mPl*zPlsqr)) \
        - c1sqr*(pow(c2,5)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + \
        30*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 3*mPl*(rt + xPl)*zPl - \
        5*c2pow4*mPl*(7*rt + xPl)*zPl - 2*c2sqr*mPl*(8*rt + xPl)*zPl + \
        2*c2cub*(2*IPlx - IPlz + 2*ITpx - ITpz - mPl*rtsqr - mTp*rtsqr + \
        4*mPl*rt*xPl - mPl*xPlsqr + mPl*yPlsqr + 2*mPl*zPlsqr) + c2*(4*IPlx - \
        3*IPlz + 4*ITpx - 3*ITpz + mPl*rtsqr + mTp*rtsqr - 6*mPl*rt*xPl - \
        3*mPl*xPlsqr + mPl*yPlsqr + 4*mPl*zPlsqr))) + \
        dc2*(2*pow(c1,6)*mPl*yPl*(2*c2*rt + zPl) + c1pow4*mPl*yPl*(-17*c2*rt \
        - 40*c2cub*rt + 4*c2*xPl + 2*zPl + 6*c2sqr*zPl) + (1 + \
        c2sqr)*mPl*yPl*(-(c2*(rt - 4*xPl)) + c2cub*(3*rt + 4*xPl) - 2*zPl \
        + 2*c2pow4*zPl) + 2*c1sqr*mPl*yPl*(10*pow(c2,5)*rt + c2*(-11*rt + \
        4*xPl) + c2cub*(-7*rt + 4*xPl) - zPl + 2*c2sqr*zPl + \
        3*c2pow4*zPl) + 2*pow(c1,5)*(IPly - IPlz + ITpy - ITpz - mPl*rtsqr - \
        mTp*rtsqr - mPl*rt*xPl - mPl*yPlsqr + c2sqr*(IPlz + ITpz + mPl*rtsqr \
        + mTp*rtsqr - 10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        2*c2*mPl*(3*rt - xPl)*zPl + mPl*zPlsqr) + c1*(2*IPly - 2*IPlz + \
        2*ITpy - 2*ITpz - mPl*rtsqr - mTp*rtsqr - 2*mPl*rt*xPl - 2*mPl*yPlsqr \
        + 2*pow(c2,6)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + \
        mPl*xPlsqr + mPl*yPlsqr) + 4*pow(c2,5)*mPl*(3*rt - xPl)*zPl - \
        8*c2cub*mPl*(2*rt + xPl)*zPl - 4*c2*mPl*(3*rt + xPl)*zPl + \
        2*mPl*zPlsqr + 2*c2pow4*(IPly + IPlz + ITpy + ITpz + mPl*rtsqr + \
        mTp*rtsqr + 13*mPl*rt*xPl + 2*mPl*xPlsqr + mPl*yPlsqr + mPl*zPlsqr) + \
        c2sqr*(4*IPly - 2*IPlz + 4*ITpy - 2*ITpz - 5*mPl*rtsqr - 5*mTp*rtsqr \
        + 12*mPl*rt*xPl + 2*mPl*xPlsqr - 2*mPl*yPlsqr + 4*mPl*zPlsqr)) + \
        c1cub*(4*IPly - 4*IPlz + 4*ITpy - 4*ITpz - 3*mPl*rtsqr - \
        3*mTp*rtsqr - 4*mPl*rt*xPl - 4*mPl*yPlsqr + 4*c2pow4*(IPlz + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) - \
        8*c2*mPl*xPl*zPl - 8*c2cub*mPl*(5*rt + xPl)*zPl + 4*mPl*zPlsqr + \
        4*c2sqr*(IPly + ITpy + mPl*(-2*rt*xPl + xPlsqr + \
        zPlsqr))))))/pow(delc,5),0,


        0,0,0,0,0,0,(4*(-(dc2*(pow(c1,7)*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr \
        - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        3*pow(c1,6)*mPl*yPl*(-4*c2*rt + zPl) + 5*c1pow4*mPl*yPl*(-(c2*rt) + \
        8*c2cub*rt + zPl + c2sqr*zPl) - (1 + c2sqr)*mPl*yPl*(-(c2*rt) + \
        c2cub*rt + zPl + 2*c2sqr*zPl + c2pow4*zPl) + \
        c1sqr*mPl*yPl*(8*c2*rt + 10*c2cub*rt - 12*pow(c2,5)*rt + zPl + \
        2*c2sqr*zPl + c2pow4*zPl) + pow(c1,5)*(2*IPly + IPlz + 2*ITpy + ITpz \
        + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + 3*mPl*xPlsqr + mPl*yPlsqr + \
        c2sqr*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + \
        mPl*xPlsqr + mPl*yPlsqr) - 16*c2*mPl*rt*zPl + 2*mPl*zPlsqr) + \
        c1cub*(4*IPly - IPlz + 4*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + \
        2*mPl*rt*xPl + 3*mPl*xPlsqr - mPl*yPlsqr - c2pow4*(IPlz + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) - \
        14*c2*mPl*rt*zPl + 40*c2cub*mPl*rt*zPl + 4*mPl*zPlsqr + \
        2*c2sqr*(2*IPly - IPlz + 2*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + \
        16*mPl*rt*xPl + mPl*xPlsqr - mPl*yPlsqr + 2*mPl*zPlsqr)) + c1*(2*IPly \
        - IPlz + 2*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + 2*mPl*rt*xPl + \
        mPl*xPlsqr - mPl*yPlsqr - pow(c2,6)*(IPlz + ITpz + mPl*rtsqr + \
        mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 2*c2*mPl*rt*zPl \
        + 10*c2cub*mPl*rt*zPl - 8*pow(c2,5)*mPl*rt*zPl + 2*mPl*zPlsqr + \
        c2pow4*(2*IPly - 3*IPlz + 2*ITpy - 3*ITpz - 3*mPl*rtsqr - 3*mTp*rtsqr \
        - 14*mPl*rt*xPl - mPl*xPlsqr - 3*mPl*yPlsqr + 2*mPl*zPlsqr) + \
        c2sqr*(4*IPly - 3*IPlz + 4*ITpy - 3*ITpz + mPl*rtsqr + mTp*rtsqr + \
        2*mPl*rt*xPl + mPl*xPlsqr - 3*mPl*yPlsqr + 4*mPl*zPlsqr)))) + \
        dc1*(2*pow(c1,6)*(c2*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr - \
        2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + mPl*(rt - xPl)*zPl) + \
        pow(c1,5)*mPl*yPl*(-((9 + 20*c2sqr)*rt) + 4*(xPl + c2*zPl)) - c1*(1 + \
        c2sqr)*mPl*yPl*((-3 - 7*c2sqr + 4*c2pow4)*rt - 4*(1 + c2sqr)*(xPl + \
        c2*zPl)) + 2*c1cub*mPl*yPl*((-3 + 13*c2sqr + 20*c2pow4)*rt + 4*(1 \
        + c2sqr)*(xPl + c2*zPl)) - (1 + c2sqr)*(-3*c2sqr*mPl*rt*zPl + \
        2*c2pow4*mPl*(rt + xPl)*zPl - mPl*(rt + 2*xPl)*zPl + c2*(-2*IPlx + \
        2*IPlz - 2*ITpx + 2*ITpz + mPl*rtsqr + mTp*rtsqr + 2*mPl*rt*xPl + \
        2*mPl*xPlsqr - 2*mPl*zPlsqr) - 2*c2cub*(IPlx - IPlz + ITpx - ITpz \
        - mPl*rtsqr - mTp*rtsqr - 3*mPl*rt*xPl - mPl*xPlsqr + mPl*zPlsqr)) + \
        c1pow4*(4*c2cub*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + \
        10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) - 6*c2sqr*mPl*(5*rt + \
        xPl)*zPl - mPl*(7*rt + 2*xPl)*zPl + 2*c2*(IPlx + IPlz + ITpx + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 7*mPl*rt*xPl + mPl*xPlsqr + 2*mPl*yPlsqr + \
        mPl*zPlsqr)) + c1sqr*(2*pow(c2,5)*(IPlz + ITpz + mPl*rtsqr + \
        mTp*rtsqr - 10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        2*c2sqr*mPl*(5*rt - 2*xPl)*zPl + 6*c2pow4*mPl*(5*rt - xPl)*zPl + \
        2*mPl*(-4*rt + xPl)*zPl + c2*(4*IPlx - 2*IPlz + 4*ITpx - 2*ITpz - \
        5*mPl*rtsqr - 5*mTp*rtsqr + 16*mPl*rt*xPl - 2*mPl*xPlsqr + \
        2*mPl*yPlsqr + 4*mPl*zPlsqr) + 4*c2cub*(IPlx + ITpx + \
        mPl*(2*rt*xPl + yPlsqr + \
        zPlsqr))))))/pow(delc,5),(-4*(dc1*(pow(c1,7)*(IPlz + ITpz + mPl*rtsqr \
        + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        3*pow(c1,6)*mPl*yPl*(-4*c2*rt + zPl) + 5*c1pow4*mPl*yPl*(-(c2*rt) + \
        8*c2cub*rt + zPl + c2sqr*zPl) - (1 + c2sqr)*mPl*yPl*(-(c2*rt) + \
        c2cub*rt + zPl + 2*c2sqr*zPl + c2pow4*zPl) + \
        c1sqr*mPl*yPl*(8*c2*rt + 10*c2cub*rt - 12*pow(c2,5)*rt + zPl + \
        2*c2sqr*zPl + c2pow4*zPl) + pow(c1,5)*(2*IPly + IPlz + 2*ITpy + ITpz \
        + mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + 3*mPl*xPlsqr + mPl*yPlsqr + \
        c2sqr*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + \
        mPl*xPlsqr + mPl*yPlsqr) - 16*c2*mPl*rt*zPl + 2*mPl*zPlsqr) + \
        c1cub*(4*IPly - IPlz + 4*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + \
        2*mPl*rt*xPl + 3*mPl*xPlsqr - mPl*yPlsqr - c2pow4*(IPlz + ITpz + \
        mPl*rtsqr + mTp*rtsqr + 30*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) - \
        14*c2*mPl*rt*zPl + 40*c2cub*mPl*rt*zPl + 4*mPl*zPlsqr + \
        2*c2sqr*(2*IPly - IPlz + 2*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + \
        16*mPl*rt*xPl + mPl*xPlsqr - mPl*yPlsqr + 2*mPl*zPlsqr)) + c1*(2*IPly \
        - IPlz + 2*ITpy - ITpz - mPl*rtsqr - mTp*rtsqr + 2*mPl*rt*xPl + \
        mPl*xPlsqr - mPl*yPlsqr - pow(c2,6)*(IPlz + ITpz + mPl*rtsqr + \
        mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 2*c2*mPl*rt*zPl \
        + 10*c2cub*mPl*rt*zPl - 8*pow(c2,5)*mPl*rt*zPl + 2*mPl*zPlsqr + \
        c2pow4*(2*IPly - 3*IPlz + 2*ITpy - 3*ITpz - 3*mPl*rtsqr - 3*mTp*rtsqr \
        - 14*mPl*rt*xPl - mPl*xPlsqr - 3*mPl*yPlsqr + 2*mPl*zPlsqr) + \
        c2sqr*(4*IPly - 3*IPlz + 4*ITpy - 3*ITpz + mPl*rtsqr + mTp*rtsqr + \
        2*mPl*rt*xPl + mPl*xPlsqr - 3*mPl*yPlsqr + 4*mPl*zPlsqr))) + \
        dc2*(4*pow(c1,7)*mPl*rt*yPl + mPl*rt*zPl - 10*c2sqr*mPl*rt*zPl + \
        5*c2pow4*mPl*rt*zPl + pow(c1,5)*mPl*yPl*((9 - 40*c2sqr)*rt + \
        4*c2*zPl) + 2*c1cub*mPl*yPl*((3 - 25*c2sqr + 10*c2pow4)*rt + \
        4*c2*(1 + c2sqr)*zPl) + c1*mPl*yPl*((1 - 10*c2sqr + 5*c2pow4)*rt + \
        4*c2*pow(1 + c2sqr,2)*zPl) + 2*pow(c1,6)*(c2*(IPlz + ITpz + mPl*rtsqr \
        + mTp*rtsqr - 10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        2*mPl*rt*zPl) + c2*(2*IPly + 2*ITpy - mPl*rtsqr - mTp*rtsqr - \
        4*mPl*rt*xPl + 2*mPl*xPlsqr + 2*mPl*zPlsqr) + c2cub*(4*IPly + \
        4*ITpy + 3*mPl*rtsqr + 3*mTp*rtsqr + 12*mPl*rt*xPl + 4*mPl*xPlsqr + \
        4*mPl*zPlsqr) + 2*pow(c2,5)*(IPly + ITpy + mPl*(xPlsqr + zPlsqr)) + \
        c1pow4*(4*c2cub*(IPlz + ITpz + mPl*rtsqr + mTp*rtsqr + \
        10*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + 9*mPl*rt*zPl - \
        40*c2sqr*mPl*rt*zPl + 2*c2*(IPly + 2*IPlz + ITpy + 2*ITpz + \
        2*mPl*rtsqr + 2*mTp*rtsqr - 22*mPl*rt*xPl + 3*mPl*xPlsqr + \
        2*mPl*yPlsqr + mPl*zPlsqr)) + c1sqr*(2*pow(c2,5)*(IPlz + ITpz + \
        mPl*rtsqr + mTp*rtsqr - 2*mPl*rt*xPl + mPl*xPlsqr + mPl*yPlsqr) + \
        6*mPl*rt*zPl - 50*c2sqr*mPl*rt*zPl + 20*c2pow4*mPl*rt*zPl + \
        4*c2cub*(IPly + IPlz + ITpy + ITpz + mPl*rtsqr + mTp*rtsqr + \
        13*mPl*rt*xPl + 2*mPl*xPlsqr + mPl*yPlsqr + mPl*zPlsqr) + c2*(4*IPly \
        + 2*IPlz + 4*ITpy + 2*ITpz + mPl*rtsqr + mTp*rtsqr - 28*mPl*rt*xPl + \
        6*mPl*xPlsqr + 2*mPl*yPlsqr + 4*mPl*zPlsqr)))))/pow(delc,5),0,


        0,0,0,0,0,0,(4*mPl*(dc1*(c1cub*yPl - 3*c1*(1 + c2sqr)*yPl + (1 + \
        c2sqr)*(c2*xPl - zPl) + 3*c1sqr*(-(c2*xPl) + zPl)) + \
        dc2*(c1cub*xPl + 3*c1sqr*c2*yPl - c2*(1 + c2sqr)*yPl + c1*(xPl - \
        3*c2sqr*xPl + 4*c2*zPl))))/pow(delc,3),(-4*mPl*(c2*(-3*dc2*xPl + \
        dc1*yPl) + c2cub*(dc2*xPl + dc1*yPl) + c1cub*(-(dc1*xPl) + \
        dc2*yPl) + dc2*zPl - 3*c2sqr*dc2*zPl + c1sqr*(-3*c2*(dc2*xPl + \
        dc1*yPl) + dc2*zPl) + c1*((1 - 3*c2sqr)*dc2*yPl + dc1*(-xPl + \
        3*c2sqr*xPl - 4*c2*zPl))))/pow(delc,3),0;
}

void InertiaModel::Gq(VectorXd &q9, VectorXd &GqVec)
{
  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  

  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double c1 = q9(6);
  double c2 = q9(7);

  double mCr = mass[0];
  double mSt = mass[1];
  double mTp = mass[2];
  double mPl = mass[3];

  double c1sqr = c1*c1;
  double c2sqr = c2*c2;

  double delc = 1+c1sqr + c2sqr;
  double delcsqr = delc*delc;

  double cth1 = cos(tht1);
  double cth2 = cos(tht2);
  double cth3 = cos(tht3); 

  double cphi1 = cos(phi1);
  double cphi2 = cos(phi2);
  double cphi3 = cos(phi3); 

  GqVec <<
          g*((lCr*mCr*cth1)/2. + lCr*mSt*cth1),g*((lCr*mCr*cth2)/2. + \
          lCr*mSt*cth2),g*((lCr*mCr*cth3)/2. + \
          lCr*mSt*cth3),(g*lSt*mSt*cphi1)/2.,(g*lSt*mSt*cphi2)/2.,(g*lSt*mSt*\
          cphi3)/2.,g*mPl*((4*c1*c2*xPl)/delcsqr - (4*c1sqr*yPl)/delcsqr + \
          (2*yPl)/delc + (2*c1*(-1 + c1sqr + c2sqr)*zPl)/delcsqr - \
          (2*c1*zPl)/delc),g*mPl*((4*c2sqr*xPl)/delcsqr - (2*xPl)/delc - \
          (4*c1*c2*yPl)/delcsqr + (2*c2*(-1 + c1sqr + c2sqr)*zPl)/delcsqr - \
          (2*c2*zPl)/delc),g*(mPl + mTp);

}

void InertiaModel::Jqtht(VectorXd &q9, MatrixXd &J_q_tht)
{

  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  
  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double c1 = q9(6);
  double c2 = q9(7);
  double z = q9(8);
  
  double c1sqr = c1*c1;
  double c2sqr = c2*c2;

  double delc = 1+c1sqr + c2sqr;

  double cth1 = cos(tht1);
  double cth2 = cos(tht2);
  double cth3 = cos(tht3); 

  double cphi1 = cos(phi1);
  double cphi2 = cos(phi2);
  double cphi3 = cos(phi3); 

  double sth1 = sin(tht1);
  double sth2 = sin(tht2);
  double sth3 = sin(tht3); 

  double sphi1 = sin(phi1);
  double sphi2 = sin(phi2);
  double sphi3 = sin(phi3); 

  MatrixXd J_eta_tht(6,3);
  MatrixXd J_eta_phi(6,6);

  J_eta_tht << 
      delc*lCr*sth1,0,0,
      -(delc*lCr*cth1),0,0,
      0,delc*lCr*sth2,0,
      0,-(delc*lCr*cth2),0,
      0,0,delc*lCr*sth3,
      0,0,-(delc*lCr*cth3) ;
  
  J_eta_phi << 
    delc*lSt*sphi1,0,0,-2*c1*b + 4*c1*rt - 2*c1*lCr*cth1 - \
    2*c1*lSt*cphi1,-2*c2*b - 4*c2*rt - 2*c2*lCr*cth1 - 2*c2*lSt*cphi1,0,

    -(delc*lSt*cphi1),0,0,2*c1*z - 2*c1*lCr*sth1 - 2*c1*lSt*sphi1,-2*rt + \
    2*c2*z - 2*c2*lCr*sth1 - 2*c2*lSt*sphi1,delc,

    0,delc*lSt*sphi2,0,-2*c1*b - 2*c1*rt - 2*sqrt3*c2*rt - \
    2*c1*lCr*cth2 - 2*c1*lSt*cphi2,-2*c2*b - 2*sqrt3*c1*rt + 2*c2*rt - \
    2*c2*lCr*cth2 - 2*c2*lSt*cphi2,0,

    0,-(delc*lSt*cphi2),0,sqrt3*rt + 2*c1*z - 2*c1*lCr*sth2 - \
    2*c1*lSt*sphi2,rt + 2*c2*z - 2*c2*lCr*sth2 - 2*c2*lSt*sphi2,delc,

    0,0,delc*lSt*sphi3,-2*c1*b - 2*c1*rt + 2*sqrt3*c2*rt - \
    2*c1*lCr*cth3 - 2*c1*lSt*cphi3,-2*c2*b + 2*sqrt3*c1*rt + 2*c2*rt - \
    2*c2*lCr*cth3 - 2*c2*lSt*cphi3,0,

    0,0,-(delc*lSt*cphi3),-(sqrt3*rt) + 2*c1*z - 2*c1*lCr*sth3 - \
    2*c1*lSt*sphi3,rt + 2*c2*z - 2*c2*lCr*sth3 - 2*c2*lSt*sphi3,delc;


  MatrixXd J_eta_phi_inv = J_eta_phi.inverse();
  MatrixXd J_phi_tht = (-J_eta_phi_inv*J_eta_tht);

  J_q_tht<<  Matrix3d::Identity(), J_phi_tht;

}