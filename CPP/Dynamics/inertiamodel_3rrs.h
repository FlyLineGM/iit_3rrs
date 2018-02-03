#ifndef INERTIAMODEL_3RRS_H
#define INERTIAMODEL_3RRS_H 1
#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Eigenvalues>

using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Infinity;
class InertiaModel
{
  double g, Eyng, density, delSt, delCr, P;
  double lCr, lSt, b, rt;
  double h0Cr, h0St, minb;
  double tTp, wTp, dt;
  double mass[4];
  double inertia[4][3];
  double xPl, yPl, zPl;
  double sqrt3;

  void Mq(VectorXd &q9, MatrixXd &MqMat);
  void Cq(VectorXd &q9, VectorXd &dq9, MatrixXd &CqMat);
  void Gq(VectorXd &q9, VectorXd &GqVec);
  void Jqtht(VectorXd &q9,  MatrixXd &J_q_tht);  

  bool rrdyad(double, double, double [2][2], double [2][2]);
  void calc();

  bool invk(double *, double *, double *);
  bool invk1(double phi[3], double omega[3], double alpha[3], double zderv[3], VectorXd &q9, VectorXd &dq9, VectorXd &ddq9);  
 public:
  InertiaModel();

  void output();
  void setlength(double *);
  void setpayload(double *, double *, double, double, double, double, double, double);
  void torquecalc(double *, double *, double *, double *, double *);
  void torquecalcmax(double *, double *, double *, double *, double &);
  void MqTerm(double *, double *, double *, double *, double, double &);
  void GqTerm(double *, double *, double *, double *, double, double &);

  double torquescan(double*, double, double*);
  void plot(double *, double *, double *, double *, double *);

  void calcMthetaEigen(VectorXd &q9, double *);
  void dynindices(double *, double *);
  bool dynScan(double *, double *);


};
#endif /*INERTIAMODEL_3RRS_H */
