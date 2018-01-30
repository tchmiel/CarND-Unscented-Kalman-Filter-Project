#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

	Tools tools;

	is_initialized_ = false;

	time_us_ = 0;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	// state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
	x_ = VectorXd(5);
	x_ << 1, 1, 1, 1, 1;

	// initial covariance matrix
	P_ = MatrixXd(5, 5);
	P_ = MatrixXd::Identity(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	// Bicycle acceleration is  0.5 to 1 m/s^2
  // https://www.fhwa.dot.gov/publications/research/safety/04103/06.cfm
	// The AASHTO Guide for the Development of Bicycle Facilities (p. 65) 
	// uses a bicycle acceleration rate of 0.5 to 1 m/sec2 
	std_a_ = 2;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.25;

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

	///* State dimension
	n_x_ = 5;

	///* Augmented state dimension
	n_aug_ = 7;

	//set measurement dimension, radar can measure r, phi, and r_dot
	n_z_ = 3;

	//LIDAR  Measurement covariance matrix
	R_laser_ = MatrixXd(2, 2);
	R_laser_ << std_laspx_ * std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	// LIDAR Measurement matrix
	H_laser_ = MatrixXd(2, 5);
	H_laser_ << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;

	//RADAR  Measurement covariance matrix
	R_radar_ = MatrixXd(n_z_, n_z_);
	R_radar_ << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	if (!is_initialized_) {
		// Initialize the state matrix, x
		x_ << 1, 1, 0, 0, 0;

		// Initialize the covariance matrix, P_
		P_ << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 0, 1, 0,
			0, 0, 0, 0, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			x_(0) = rho * cos(phi);
			x_(1) = rho * sin(phi);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
			/**
			Initialize state.
			*/
			x_(0) = meas_package.raw_measurements_[0];
			x_(1) = meas_package.raw_measurements_[1];
		}

		time_us_ = meas_package.timestamp_;

		//create augmented sigma points 
		Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
		Xsig_aug_.fill(0.0);

		//create matrix with predicted sigma points as columns
		Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
		Xsig_pred_.fill(0.0);

		//create vector for weights
		weights_ = VectorXd(2 * n_aug_ + 1);
		weights_.fill(0.0);

		///* mean predicted measurement
		z_pred_ = VectorXd(n_z_);
		z_pred_.fill(0.0);

		///* measurement covariance matrix S
		S_ = MatrixXd(n_z_, n_z_);
		S_.fill(0.0);

		//create matrix for sigma points in measurement space
		Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
		Zsig_.fill(0.0);

		// spreading parameter
		lambda_ = 3 - n_aug_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}

	//dt - expressed in seconds
	float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	Prediction(dt);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR  && use_radar_) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER  && use_laser_) {
		UpdateLidar(meas_package);
	}

	// print the output
	//cout << "x_ = " << x_ << endl;
	//cout << "P_ = " << P_ << endl;
}


void UKF::PredictMeanAndCovariance() {

	// set weights
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i = 1; i < 2 * n_aug_ + 1; i++) {
		//2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda_);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//iterate over sigma points
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		
		//angle normalization
		x_diff(3) = tools.NormalizeAngle(x_diff(3));

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}

}

void UKF::SigmaPointPrediction(double delta_t) {

	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug_(0, i);
		double p_y = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double yaw = Xsig_aug_(3, i);
		double yawd = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
		}
		else {
			px_p = p_x + v * delta_t*cos(yaw);
			py_p = p_y + v * delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd * delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	//print result
	//std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;
}

void UKF::AugmentedSigmaPoints() {

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	// Create Augmented Sigma Points
	//    Lession 7, section 18:  Augmentation Assignment 2
	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	float spread_factor = sqrt(lambda_ + n_aug_);
	Xsig_aug_.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug_.col(i + 1) = x_aug + spread_factor * L.col(i);
		Xsig_aug_.col(i + 1 + n_aug_) = x_aug - spread_factor * L.col(i);
	}
	
	//print result
	//std::cout << "Xsig_aug_ = " << std::endl << Xsig_aug_ << std::endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	VectorXd z = meas_package.raw_measurements_;

	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;

	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;

	// TODO:  Calculate the Lidar NIS
}

void UKF::UpdateState(MeasurementPackage meas_package) {
	//calculate cross correlation matrix
	MatrixXd Tc = MatrixXd(n_x_, n_z_);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
		//2n+1 sigma points

		//residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;
		
		//angle normalization
		z_diff(1) = tools.NormalizeAngle(z_diff(1));
		
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		
		//angle normalization
		x_diff(3) = tools.NormalizeAngle(x_diff(3));
		
		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S_.inverse();

	//residual
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff = z - z_pred_;

	//angle normalization
	z_diff(1) = tools.NormalizeAngle(z_diff(1));

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S_ * K.transpose();

	//print result
	//std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	//std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
}

void UKF::PredictRadarMeasurement() {
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		//2n+1 sigma points

		// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig_(0, i) = sqrt(p_x*p_x + p_y * p_y);							//r
		Zsig_(1, i) = atan2(p_y, p_x);										//phi
		Zsig_(2, i) = (p_x*v1 + p_y * v2) / sqrt(p_x*p_x + p_y * p_y);		//r_dot
	}

	//mean predicted measurement
	z_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
	}

	//innovation covariance matrix S
	S_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd z_diff = Zsig_.col(i) - z_pred_;

		//angle normalization
		z_diff(1) = tools.NormalizeAngle(z_diff(1));

		S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S_ = S_ + R_radar_;

	//print result
	//std::cout << "z_pred: " << std::endl << z_pred_ << std::endl;
	//std::cout << "S: " << std::endl << S_ << std::endl;
}
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	PredictRadarMeasurement();
	UpdateState(meas_package);
	
	//calculate the radar NIS
}
