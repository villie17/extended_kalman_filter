#include "kalman_filter.h"
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#define MY_DEBUG
#undef MY_DEBUG

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long  x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  double phi =  sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
  double rho = atan2(x_[1], x_[0]);
  double phi_dot = (x_[0]*x_[2] + x_[1]*x_[3])/phi;

  VectorXd hx(3);
  hx << phi, rho, phi_dot;


#ifdef MY_DEBUG
  for (int i=0; i<3; i++){
	  cout << std::setprecision(5) << z[i] << std::setprecision(5)<< "\t:" << hx[i]<<'|';
  }
  cout << endl;
#endif

  VectorXd y = z - hx;
  while (y[1] >= M_PI ){
  		  y[1]-= 2*M_PI;
  }
  while (y[1] <= -M_PI ){
	  y[1]+= 2*M_PI;
  }
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long  x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
