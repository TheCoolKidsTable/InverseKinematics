/*
 * RBDL - Rigid Body Dynamics Library
 * Copyright (c) 2011-2016 Martin Felis <martin@fysx.org>
 *
 * Licensed under the zlib license. See LICENSE for more details.
 */

#include <iostream>
#include <chrono>
#include <thread>

#include <rbdl/rbdl.h>
#include <rbdl/rbdl_utils.h>


#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include <numeric>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include <math.h>
#include <nlopt.hpp>
// #include <nlopt.h>

#include "rotation.hpp"

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

typedef struct {
	Model* model;
	Eigen::Affine3d target;
	int body_id;
} params;


void forwardKinematics(Math::VectorNd & q, Math::Vector3d & xe, Math::Matrix3d & R, Model * model, int body_id) {
    Math::Vector3d base;
    base << 0,0,0;
    // position
	xe = CalcBaseToBodyCoordinates(*model, q, body_id, base);
    // orientation
    R = CalcBodyWorldOrientation(*model, q, body_id);
}

double cost(const std::vector<double> & q_nlopt, std::vector<double> &grad, void *data)
{
	params* p = reinterpret_cast<params*>(data);
	Model* model = p->model;
	int body_id = p->body_id;

	Math::VectorNd q;
	q.resize(q_nlopt.size());
	for(int i = 0; i < q_nlopt.size(); i++) {
		q(i) = q_nlopt[i];
	}
	Math::Vector3d xe;
	Math::Matrix3d R;
	Math::Matrix3d Rd;

	//  calculate current end effector position
	forwardKinematics(q, xe, R, model, body_id);

	// decompose the euler angles from current forward kinematics
	Eigen::MatrixXd Theta_e = R.eulerAngles(2,1,0);
	// decompose the euler angles from target orientation
	Eigen::MatrixXd Theta_d = p->target.linear().eulerAngles(2,1,0);

	// Tuning
	Eigen::MatrixXd K(3,3);
	K.diagonal() << 0.00001, 0.00001, 0.00001;

	// position
	Eigen::Vector3d xd = p->target.translation();
	Eigen::MatrixXd position_cost = xe - xd;
	// orientation
	Eigen::MatrixXd orientation_cost = (Theta_d - Theta_e).transpose()*K*(Theta_d - Theta_e);
	double cost =  position_cost.norm() + orientation_cost(0,0);
	return cost;
}

int main (int argc, char* argv[]) {
	rbdl_check_api_version (RBDL_API_VERSION);

	Model* model = new Model();

	if (argc != 2) {
		std::cerr << "Error: not right number of arguments." << std::endl;
		std::cerr << "usage: " << argv[0] << " <model.urdf>" << std::endl;
		exit(-1);
	}

	if (!Addons::URDFReadFromFile (argv[1], model, false)) {
		std::cerr << "Error loading model " << argv[1] << std::endl;
		abort();
	}

	std::cout << "Degree of freedom overview:" << std::endl;
	std::cout << Utils::GetModelDOFOverview(*model);

	std::cout << "Model Hierarchy:" << std::endl;
	std::cout << Utils::GetModelHierarchy(*model);

	VectorNd Q = VectorNd::Zero (model->q_size);
	VectorNd QDot = VectorNd::Zero (model->qdot_size);
	VectorNd Tau = VectorNd::Zero (model->qdot_size);
	VectorNd QDDot = VectorNd::Zero (model->qdot_size);

	std::cout << "Forward Dynamics with q, qdot, tau set to zero:" << std::endl;
	ForwardDynamics (*model, Q, QDot, Tau, QDDot);

	std::cout << QDDot.transpose() << std::endl;

	std::cout << "************************ Forward Kinematics ************************" << std::endl;

	Math::Vector3d base;
	base << 0,0,0;
	Math::Vector3d foot_pos;
	foot_pos = CalcBaseToBodyCoordinates(*model, Q, 16, base);
	Math::Matrix3d R;
	R = CalcBodyWorldOrientation(*model, Q, 16);

	std::cout << "Forward Kinematics result pos: "<< foot_pos << std::endl;
	std::cout << "Forward Kinematics result orientation: "<< R << std::endl;
	Eigen::Matrix<double, Eigen::Dynamic, 1> Theta;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Rn = R;
	rot2rpy(Rn,Theta);
	std::cout << "Forward Kinematics result euler angles: "<< Theta << std::endl;

	std::cout << "************************ Inverse Kinematics ************************" << std::endl;
	Math::VectorNd Qres;
	std::vector<unsigned int> body_ids;
	body_ids.push_back(5);
	std::vector<Math::Vector3d> body_points;
	Math::Vector3d body_point;
	body_point << 0,0,0;
	body_points.push_back(body_point);
	std::vector<Math::Vector3d> target_positions;
	Math::Vector3d target_position;
	target_position << 0,0,0;
	target_positions.push_back(target_position);

	auto t_start = std::chrono::high_resolution_clock::now();
	InverseKinematics(*model,Q,body_ids,body_points,target_positions,Qres);
	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	std::cout << "Time taken for IK [ms] : " << elapsed_time_ms  << std::endl;
	std::cout << "Inverse Kinematics result: "<< Qres << std::endl;
	foot_pos = CalcBaseToBodyCoordinates(*model, Qres, 5, base);
	std::cout << "Forward Kinematics result for IK solution: "<< foot_pos << std::endl;

	std::cout << "************************ Inverse Kinematics with Constraint Set ************************" << std::endl;

	InverseKinematicsConstraintSet CS;

	std::vector<Math::Matrix3d> target_orientations;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  target_orientation;

	roty(M_1_PI/2,target_orientation);

	target_orientation = target_orientation*Rn;

	rot2rpy(target_orientation,Theta);
	std::cout << "Forward Kinematics result euler angles: "<< Theta << std::endl;

	target_orientations.push_back(target_orientation);

	int body_id = 4;

	std::cout << "Target foot pos: "<< foot_pos << std::endl;
	// CS.AddOrientationConstraint(body_id,target_orientation);
	// CS.AddFullConstraint(body_id,body_point,target_position,R);
	CS.AddPointConstraint(body_id,body_point,target_position);
	CS.constraint_tol = 1e-12;
	CS.max_steps= 300;
	CS.lambda = 1e-9;

	auto t_start2 = std::chrono::high_resolution_clock::now();
	InverseKinematics(*model, Q, CS, Qres);
	auto t_end2 = std::chrono::high_resolution_clock::now();
	auto elapsed_time_ms2 = std::chrono::duration<double, std::milli>(t_end2-t_start2).count();
	std::cout << "Time taken for IK with constraints [ms] : " << elapsed_time_ms2  << std::endl;
	std::cout << "Inverse Kinematics result: "<< Qres << std::endl;
	foot_pos = CalcBaseToBodyCoordinates(*model, Qres, 5, base);
	R = CalcBodyWorldOrientation(*model, Qres, 5);

	std::cout << "Forward Kinematics position result for IK solution: "<< foot_pos << std::endl;
	std::cout << "Forward Kinematics orientation result for IK solution: "<< R << std::endl;
	Math::Vector3d xe;

	std::cout << "************************ NLopt Forward Kinematics ************************" << std::endl;

	Q.setZero();
	forwardKinematics(Q, xe, R, model, 5);

	std::cout << "Forward Kinematics position result : " << std::endl << xe << std::endl;
	std::cout << "Forward Kinematics orientation result : " << std::endl<< R << std::endl;
	std::cout << "Forward Kinematics theta result : "<<std::endl<< R.eulerAngles(2,1,0) << std::endl;

	Eigen::Affine3d target;
	Eigen::Vector3d orientation_target;
	Eigen::Vector3d position_target;
	position_target <<  0.0550006, 5.6075e-05, -0.4;
	orientation_target << 0,1.571,0;
	target.translation() = position_target;
	target.linear() = R;

	params data[3] = { model, target, 5};

	nlopt::opt opt (nlopt::LN_COBYLA, model->q_size);

	opt.set_min_objective(cost, &data[0]);
	opt.set_xtol_rel(1e-5);
	std::vector<double> q_nlopt(model->q_size);
	// Initial conditions
	std::fill(q_nlopt.begin(), q_nlopt.end(), 0);
	double minf;
	// Run the optimizer
	auto t_start_3 = std::chrono::high_resolution_clock::now();
	std::cout << opt.optimize(q_nlopt, minf) << std::endl;
	auto t_end_3 = std::chrono::high_resolution_clock::now();
	auto elapsed_time_ms_3 = std::chrono::duration<double, std::milli>(t_end_3-t_start_3).count();
	std::cout << "Time taken for IK nlopt [ms] : " << elapsed_time_ms_3  << std::endl;

	// Print the results
	for(int i = 0; i < q_nlopt.size(); i++) {
		std::cout << "q indx : " << i << " : " << q_nlopt[i] << std::endl;
	}

	Math::VectorNd q;
	q.resize(q_nlopt.size());
	for(int i = 0; i < q_nlopt.size(); i++) {
		q(i) = q_nlopt[i];
	}

	xe = CalcBaseToBodyCoordinates(*model, q, 5, base);
    // orientation
    R = CalcBodyWorldOrientation(*model, q, 5);

	std::cout << "Inverse Kinematics position result : "<< std::endl<< xe << std::endl;
	std::cout << "Inverse Kinematics orientation result : "<<std::endl<< R << std::endl;
	std::cout << "Inverse Kinematics theta result : "<<std::endl<< R.eulerAngles(2,1,0) << std::endl;


	delete model;
	return 0;
}



