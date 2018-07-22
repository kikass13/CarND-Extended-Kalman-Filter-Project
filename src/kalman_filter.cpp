#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in)
{
	this->x_ = x_in;
	this->P_ = P_in;
	this->F_ = F_in;
	//these come from outside before update!
	//this->H_ = H_in;
	//this->R_ = R_in;
	//this->Q_ = Q_in;
	//set our n dims
	this->n = this->x_.size();
}

void KalmanFilter::Predict() 
{
	/**
	TODO:
		* predict the state
	*/
	this->x_ = this->F_ * this->x_;
	this->P_ = this->F_ * this->P_ * this->F_.transpose() + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
	/**
	TODO:
		* update the state by using Kalman Filter equations
	*/
	VectorXd z_pred = this->H_ * this->x_;
	VectorXd y = z - z_pred;
	//kalman gain
	MatrixXd S = this->H_ * this->P_ * this->H_.transpose() + this->R_;
	MatrixXd K = this->P_ * this->H_.transpose() * S.inverse();
	//posterior state
	this->x_ = this->x_ + (K * y);
	MatrixXd I = MatrixXd::Identity(this->n, this->n);
	//update posterior cov via kalman gain and previous uncertainty
	this->P_ = (I - K * this->H_) * this->P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
	/**
	TODO:
		* update the state by using Extended Kalman Filter equations
	*/
	//grab important measurement values
	float rho = sqrt(this->x_(0) * this->x_(0) + this->x_(1) * this->x_(1));
	float phi = atan2(this->x_(1), this->x_(0));
	float rho_dot = (this->x_(0) * this->x_(2) + this->x_(1) * this->x_(3)) / rho;
	VectorXd z_pred(3);
	z_pred << rho, phi, rho_dot;
	VectorXd y = z - z_pred;

	// Normalize angles
	if(y(1) < -M_PI)	y(1) += 2 * M_PI;
	if(y(1) > M_PI)		y(1) -= 2 * M_PI;

	//kalman gain
	MatrixXd S = this->H_ * this->P_ * this->H_.transpose() + this->R_;
	MatrixXd K = this->P_ * this->H_.transpose() * S.inverse();
	//posterior state
	this->x_ = this->x_ + (K * y);
	std::cout << this->n << "\n";
	MatrixXd I = MatrixXd::Identity(this->n, this->n);
	//update posterior cov via kalman gain and previous uncertainty
	this->P_ = (I - K * this->H_) * this->P_;
}
