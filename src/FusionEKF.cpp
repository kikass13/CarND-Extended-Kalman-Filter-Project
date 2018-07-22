#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() 
{
	this->is_initialized_ = false;

	this->previous_timestamp_ = 0;

	// initializing matrices
	this->R_laser_ = MatrixXd(2, 2);
	this->R_radar_ = MatrixXd(3, 3);
	this->H_laser_ = MatrixXd(2, 4);
	this->Hj_ = MatrixXd(3, 4);

	//measurement covariance matrix - laser
	this->R_laser_ << 0.0225, 0,
				0, 0.0225;

	//measurement covariance matrix - radar
	this->R_radar_ << 0.09, 0, 0,
				0, 0.0009, 0,
				0, 0, 0.09;

	/**
	TODO:
		* Finish initializing the FusionEKF.
		* Set the process and measurement noises
	*/
	this->Hj_ << 	1, 1, 0, 0,
		 	1, 1, 0, 0,
			1, 1, 1, 1;

	this->H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (! this->is_initialized_) 
	{
		std::cout << "Init EKF...\n";
		/**
		TODO:
			* Initialize the state ekf_.x_ with the first measurement.
			* Create the covariance matrix.
			* Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
		// first measurement
		cout << "EKF: " << endl;
		Eigen::VectorXd x = VectorXd(4);
		x << 1, 1, 1, 1;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
		{
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = measurement_pack.raw_measurements_(0);
			float phi = measurement_pack.raw_measurements_(1);
			float rho_dot = measurement_pack.raw_measurements_(2);
			x(0) = rho * cos(phi);
			x(1) = rho * sin(phi);      
			x(2) = rho_dot * cos(phi);
			x(3) = rho_dot * sin(phi);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
		{
			/**
			Initialize state.
			*/
			x(0) = measurement_pack.raw_measurements_(0);
			x(1) = measurement_pack.raw_measurements_(1);
		}

		//state covariance matrix P
		MatrixXd P = MatrixXd(4, 4);
		P << 	1, 0, 0, 0,
				 		0, 1, 0, 0,
						0, 0, 1000, 0,
						0, 0, 0, 1000;

		//the initial transition matrix F_
		MatrixXd F = MatrixXd(4, 4);
		F << 	1, 0, 1, 0,
				 		0, 1, 0, 1,
						0, 0, 1, 0,
						0, 0, 0, 1;

		//init our ekf
		this->ekf_.Init(x, P, F);


		// done initializing, no need to predict or update
		this->is_initialized_ = true;
		
		// Store the timestamp of the first measurement
		this->previous_timestamp_ = measurement_pack.timestamp_;
		std::cout << "Init Done ...\n";
		return;
	}

	//compute the time dt[s] elapsed
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	this->previous_timestamp_ = measurement_pack.timestamp_;


	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	//shortages for different pow of dt
	float dt2 = dt * dt;
	float dt3 = dt2 * dt;
	float dt4 = dt2 * dt2;

	//Modify the F matrix so that the time is integrated
	this->ekf_.F_(0, 2) = dt;
	this->ekf_.F_(1, 3) = dt;

	//set the acceleration noise components
	float noise_ax = 9;
	float noise_ay = 9;

	//set the process covariance matrix Q
	this->ekf_.Q_ = MatrixXd(4, 4);
	this->ekf_.Q_ << 	dt4/4 * noise_ax, 0, dt3/2 * noise_ax, 0,
						0, dt4/4 * noise_ay, 0, dt3/2 * noise_ay,
						dt3/2 * noise_ax, 0, dt2 * noise_ax, 0,
						0, dt3/2*noise_ay, 0, dt2*noise_ay;

	/**
	 TODO:
		 * Update the state transition matrix F according to the new elapsed time.
			- Time is measured in seconds.
		 * Update the process noise covariance matrix.
		 * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	 */

	//do ekf prediction step
	this->ekf_.Predict();


	/*****************************************************************************
	 *  Update
	 ****************************************************************************/

	/**
	 TODO:
		 * Use the sensor type to perform the update step.
		 * Update the state and covariance matrices.
	 */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
	{
		// Radar updates
		this->ekf_.H_ = tools.CalculateJacobian(this->ekf_.x_);
		this->ekf_.R_ = this->R_radar_;
		this->ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} 
	else 
	{
		// Laser updates
		this->ekf_.H_ = this->H_laser_;
		this->ekf_.R_ = this->R_laser_;
		this->ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << this->ekf_.x_ << endl;
	cout << "P_ = " << this->ekf_.P_ << endl;
}
