#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth){
  /**
	TODO:
   * Calculate the RMSE here.
   */


  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // TODO: YOUR CODE HERE

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here

  if (estimations.size() == 0
	   || estimations.size() != ground_truth.size() ){
	  return rmse;
  }

  //accumulate squared residuals

  for(unsigned int i=0; i < estimations.size(); ++i){
	  // ... your code here
	  VectorXd diff = (estimations[i] - ground_truth[i]);
	  VectorXd diff_sqrd = (diff.array() * diff.array());
	  rmse += diff_sqrd;
  }

  //calculate the mean
  // ... your code here
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();

  //calculate the squared root
  // ... your code here

  //return the result
  return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
	* Calculate a Jacobian here.
  */

  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //TODO: YOUR CODE HERE
  double px2_py2 = px*px + py*py;
  double px2_py2_sqr = sqrt(px2_py2);
  double px2_py2_3by2 = px2_py2_sqr*px2_py2;
  double vxpy_vypx = vx*py - vy * px;
  double vypx_vxpy = -vx*py + vy * px;
  double px_by_px2_py2_sqr = px/px2_py2_sqr;
  double py_by_px2_py2_sqr = py/px2_py2_sqr;


  if (px2_py2 < 0.0001){
	  std::cerr << "Div by zero" << std::endl;
	  return Hj;
  }

  Hj << px_by_px2_py2_sqr, py_by_px2_py2_sqr, 0 , 0
		  ,-py/px2_py2, px/px2_py2, 0, 0
		  ,py*(vxpy_vypx)/px2_py2_3by2, px*(vypx_vxpy)/px2_py2_3by2, px_by_px2_py2_sqr, py_by_px2_py2_sqr;



  //check division by zero

  //compute the Jacobian matrix

  return Hj;
}
