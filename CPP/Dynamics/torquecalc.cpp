// Function to return actuator torques

#include "inertiamodel_3rrs.h"
void InertiaModel::torquecalcmax(double phi[3], double omega[3], double alpha[3], double zderv[3], double &torque)
{

  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);
  VectorXd GqVec(9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);
  VectorXd Mqterm(3);
  VectorXd Gterm(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);


  /********Torque calculation********/

  invk1(phi, omega, alpha, zderv,q9, dq9, ddq9);
  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);
  Gq(q9, GqVec);
  Jqtht(q9, J_q_tht);
  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
  Toutheta = J_q_tht_trans*MqMat*J_q_tht*ddtht + J_q_tht_trans*CqMat*J_q_tht*dtht + J_q_tht_trans*GqVec;



//  std::cout<<Gterm;

  torque = Toutheta.lpNorm<Infinity>();
}

void InertiaModel::torquecalc(double phi[3], double omega[3], double alpha[3], double zderv[3], double torque[3])
{

  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);
  VectorXd GqVec(9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);
  VectorXd Mqterm(3);
  VectorXd Gterm(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);


  /********Torque calculation********/

  invk1(phi, omega, alpha, zderv,q9, dq9, ddq9);
  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);
  Gq(q9, GqVec);
  Jqtht(q9, J_q_tht);
  
  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
  Toutheta = J_q_tht_trans*MqMat*J_q_tht*ddtht + J_q_tht_trans*CqMat*J_q_tht*dtht +  J_q_tht_trans*GqVec;

  Mqterm = J_q_tht_trans*MqMat*J_q_tht*ddtht;
  Gterm = J_q_tht_trans*GqVec;

  for(int i = 0; i<3; i++)
  {
    torque[i] = Toutheta(i);
  }

}


void InertiaModel::GqTerm(double phi[3], double omega[3], double alpha[3], double zderv[3], double leg, double &torque)
{

  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);
  VectorXd GqVec(9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);
  VectorXd Mqterm(3);
  VectorXd Gterm(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);


  /********Torque calculation********/

  invk1(phi, omega, alpha, zderv,q9, dq9, ddq9);
  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);
  Gq(q9, GqVec);
  Jqtht(q9, J_q_tht);
  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
//  Toutheta = J_q_tht_trans*MqMat*J_q_tht*ddtht + J_q_tht_trans*CqMat*J_q_tht*dtht +  J_q_tht_trans*GqVec;

//  Mqterm = J_q_tht_trans*MqMat*J_q_tht*ddtht;
  Gterm = J_q_tht_trans*GqVec;

  torque = Gterm(leg);
    
}

void InertiaModel::MqTerm(double phi[3], double omega[3], double alpha[3], double zderv[3], double leg, double &torque)
{

  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);
  VectorXd GqVec(9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);
  VectorXd Mqterm(3);
  VectorXd Gterm(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);


  /********Torque calculation********/

  invk1(phi, omega, alpha, zderv,q9, dq9, ddq9);
  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);
  Gq(q9, GqVec);
  Jqtht(q9, J_q_tht);
  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
//The contribution from mass/inertia term for torque
  Mqterm = J_q_tht_trans*MqMat*J_q_tht*ddtht;

  torque = Mqterm(leg);
    
}
