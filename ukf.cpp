#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 30;
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // std_yawdd_ = 30;
  std_yawdd_ = 3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  // initial state vector
  // px, py, vel, yaw, yawd
  x_ = VectorXd(n_x_);

  //define spreading parameter
  // ????????? Should this be 3 - n_x ????
  lambda_ = 3 - n_aug_;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // set weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  cout << "weights_(0)" << weights_(0) << endl;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights_
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
    cout << "weights_(i)" << weights_(i) << endl;
  }

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  // Doesn't seem to need initialization, but just to be safe..
  Xsig_pred_.fill(0.0);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  double previous_timestamp;

  // Code from my EKF Project
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ukf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */


    // !!!! WALKTHROUGH SUGGESTS THAT state x and covariance P_ get initialized here !!!!


    // first measurement
    // cout << "EKF: " << endl;
    // I changed these to zero for ukf. Not sure why 1 would be better?
    x_ << 0, 0, 0, 0, 0;


    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // cout << "RADAR Measurement (Initial)" << endl;

      // ??????? Need to intialize differently for Radar and Lidar ????
      // Need to use manufacturer information
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho_measured = measurement_pack.raw_measurements_[0];
      double phi_measured = measurement_pack.raw_measurements_[1];
      double rhodot_measured = measurement_pack.raw_measurements_[2];

      // ?????? REVIEW IF ANYTHING NEEDS TO BE DONE FOR rhodot_measured
      // From radar angle is from your perpective, from model it is from objects perpspective
      x_ << rho_measured*cos(phi_measured), rho_measured*sin(phi_measured), rhodot_measured, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // cout << "LASER Measurement (Initial)" << endl;

      // ??????? Need to intialize differently for Radar and Lidar ????
      // Need to use manufacturer information
      P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      //set the state with the initial location and zero velocity
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;

    }

    previous_timestamp = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
  */
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp) / 1000000.0;   //dt - expressed in seconds
  previous_timestamp = measurement_pack.timestamp_;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated (section 8 of lesson 5)
  // ukf_.F_(0, 2) = dt;
  // ukf_.F_(1, 3) = dt;

  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, dt, 0,
       0, 1, 0, dt,
       0, 0, 1, 0,
       0, 0, 0, 1;

  //set the process covariance matrix Q (section 9 of lesson 5)
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  MatrixXd Q = MatrixXd(4, 4);
  Q <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
        0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
        0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  cout << "Calling Prediction()" << endl;
  Prediction(dt);


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
    cout << "Calling UpdateRadar()" << endl;
    UKF::UpdateRadar(measurement_pack);
    cout << "Returned from UpdateRadar()" << endl;
  } else {
    // Laser updates
    cout << "Calling UpdateLidar()" << endl;
    UKF::UpdateLidar(measurement_pack);
    cout << "Returned from UpdateLidar()" << endl;
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // From Lesson 7, section 18 Augmentation Assignment 2

   // CREATE AUGMENTED SIGMA POINTS

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/
 
  //create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  MatrixXd P_aug;
  P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  cout << "create augmented sigma points" << endl;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  /*******************************************************************************
   * Student part end
   ******************************************************************************/



  // From Lesson 7, section 21 Sigma Point Prediction Assignment 2
  // PREDICT SIGMA POINTS (Predict Sigma Points)

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/

  //predict sigma points
  cout << "poredict sigma points" << endl;
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into correct column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /*******************************************************************************
   * Student part end
   ******************************************************************************/

  // From Lesson 7, section 24 Predicted Mean and Variance Assignment 2
  // PREDICT MEAN & VARIANCE (Predict Mean and Variance)

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/
  
    //predicted state mean
    cout << "predict state mean" << endl;
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }
  
    //predicted state covariance matrix
    cout << "predict state covariance" << endl;
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
  
      cout << "i:" << i << endl;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      float factor = x_diff(3)/ M_PI;
      if (factor > 1)
        cout << "factor:" << factor << endl;
        x_diff(3) = x_diff(3) - (factor * M_PI);
      if (factor < 1)
        cout << "factor:" << factor << endl;
        x_diff(3) = x_diff(3) + (factor * M_PI);

/****
      while (x_diff(3)> M_PI) {
        cout << "x_diff(3)" << x_diff(3) << endl;
        cout << "Xsig_pred_.col(i):" << Xsig_pred_.col(i) << endl;
        cout << "x_:" << x_ << endl;
        x_diff(3)-=2.*M_PI;
      }
      while (x_diff(3)<-M_PI) {
        cout << "x_diff(3)" << x_diff(3) << endl;
        cout << "Xsig_pred_.col(i):" << Xsig_pred_.col(i) << endl;
        cout << "x_:" << x_ << endl;
        x_diff(3)+=2.*M_PI;
      }
****/
  
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
      
    }
  
  
  /*******************************************************************************
   * Student part end
   ******************************************************************************/
    cout << "End of Prediction()" << endl;
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.


  !!!!!! NEED TO CALCULATE lidar NIS !!!!


  */

  //set measurement dimension, lidar can measure ?????
  int n_z = 2;

  VectorXd z = VectorXd(n_z);
  z = measurement_pack.raw_measurements_;

  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // Expanded to 5 wide with 0's for state space of 5
  MatrixXd H_laser = MatrixXd(2, 5);
  H_laser << 1, 0, 0, 0, 0,
             0, 1, 0, 0, 0;

  //measurement covariance matrix - laser
  // ????? R is 1X1 for KF and 3X3 for radar for EKF ????
  MatrixXd R_laser = MatrixXd(2, 2);
  R_laser << std_laspx_*std_laspx_, 0,
             0, std_laspy_*std_laspy_;

  MatrixXd R_radar = MatrixXd(3, 3);
  R_radar << std_radr_*std_radr_, 0, 0,
             0, std_radphi_*std_radphi_, 0,
             0, 0, std_radrd_*std_radrd_;


  // From section 7 of lesson 5
  VectorXd y = z - H_laser * x_;
  MatrixXd Ht = H_laser.transpose();
  MatrixXd S = H_laser * P_ * Ht + R_laser;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_laser) * P_;
}




/*
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // From Lesson 7, section 27 Predict Radar Measurement Assignment 2
  // PREDICT RADAR MEASUREMENTS (Predict Radar Measurement)

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  VectorXd z = VectorXd(n_z);
  z = measurement_pack.raw_measurements_;

  //transform sigma points into measurement space
  cout << "HERE A" << endl;
  VectorXd Zsig = VectorXd(n_z, 2*n_aug_+1);
  cout << "HERE B" << endl;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  /*******************************************************************************
   * Student part end
   ******************************************************************************/

  // From Lesson 7, section 30 UKF Update Assignment 2
  // UKF UPDATE (UKF Update)

  /*******************************************************************************
   * Student part begin
   ******************************************************************************/

  //calculate cross correlation matrix
  VectorXd Tc = VectorXd(n_z, 2*n_aug_+1);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /*******************************************************************************
   * Student part end
   ******************************************************************************/
}

