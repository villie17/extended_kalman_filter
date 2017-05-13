#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

#define MY_DEBUG
#undef MY_DEBUG

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
		  	  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
		  	  0, 1, 0, 1,
  			  0, 0, 1, 0,
  			  0, 0, 0, 1;

  //set the acceleration noise components
  H_laser_ << 1, 0, 0, 0,
  	0, 1, 0, 0;



}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    	float x = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
    	float y = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
    	ekf_.x_ << x, y, 0, 0;
    	previous_timestamp_ = measurement_pack.timestamp_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    	previous_timestamp_ = measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = ekf_.F_(1,3) = dt;

  //2. Set the process covariance matrix Q
  double dt_4_by_4 = pow(dt, 4)/4;
  double dt_3_by_2 = pow(dt, 3)/2;

  double noise_ax = 9;
  double noise_ay = 9;


  ekf_.Q_ = MatrixXd(4,4);

  ekf_.Q_ << noise_ax *  dt_4_by_4 , 0 , noise_ax * dt_3_by_2 , 0
		  , 0 , noise_ay * dt_4_by_4 , 0 , noise_ay * dt_3_by_2
		  , noise_ax * dt_3_by_2 , 0 , noise_ax * dt * dt , 0
		  , 0 , noise_ay * dt_3_by_2 , 0 , noise_ay * dt * dt;


  ekf_.Predict();
#ifdef MY_DEBUG
  cout << "Predicted:" << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "__"<< endl;
#endif

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
	  Tools t;
	  ekf_.H_ = t.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  VectorXd z(3);
	  z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1] , measurement_pack.raw_measurements_[2];
	  ekf_.UpdateEKF(z);

  } else {
	  // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  VectorXd z(2);
	  z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
	  ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
