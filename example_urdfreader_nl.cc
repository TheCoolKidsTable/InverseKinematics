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

#include <nlopt.hpp>
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
#include <nlopt.hpp>
#include <vector>

#include "rotation.hpp"

#ifndef RBDL_BUILD_ADDON_URDFREADER
#error "Error: RBDL addon URDFReader not enabled."
#endif

#include <rbdl/addons/urdfreader/urdfreader.h>

using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

// typedef struct {
// 	Model* model;
// 	Eigen::Vector3d target;
// 	int body_id;
// } Data;

// q - joints
// xe - pos of end effector
//  R - orientation of end effector
// void forwardKinematics(Math::VectorNd & q, Math::Vector3d & xe, Math::Matrix3d & R, Model * model, int body_id) {
//     Math::Vector3d base;
//     base << 0,0,0;
//     // position
// 	std::cout << "here" << std::endl;
// 	xe = CalcBaseToBodyCoordinates(*model, q, body_id, base);
//     // orientation
//     R = CalcBodyWorldOrientation(*model, q, body_id);
// }

// double cost(const std::vector<double> & q_nlopt, std::vector<double> &grad, void *data)
// {
// 	Data* d = reinterpret_cast<Data*>(data);
// 	Model* model = d->model;
// 	int body_id = d->body_id;

// 	Math::VectorNd q(q_nlopt.data());
// 	Math::Vector3d xe;
// 	Math::Matrix3d R;

// 	// calculate current end effector position
// 	forwardKinematics(q, xe, R, model, body_id);

// 	// x desired
// 	Eigen::Vector3d xd = d->target;
// 	Eigen::Vector3d cost = xe - xd;
// 	return cost.norm();
// }

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

	std::cout << "Forward Kinematics" << std::endl;

	Math::Vector3d base;
	base << 0,0,0;

	Math::Vector3d foot_pos;

	foot_pos = CalcBaseToBodyCoordinates(*model, Q, 5, base);

	std::cout << "Forward Kinematics result: "<< foot_pos << std::endl;

	std::cout << "Inverse Kinematics" << std::endl;
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

	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	auto t_start = std::chrono::high_resolution_clock::now();
	InverseKinematics(*model,Q,body_ids,body_points,target_positions,Qres);
	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	std::cout << "Time taken for IK [ms] : " << elapsed_time_ms  << std::endl;

	std::cout << "Inverse Kinematics result: "<< Qres << std::endl;

	foot_pos = CalcBaseToBodyCoordinates(*model, Qres, 5,base);
	std::cout << "Forward Kinematics result for IK solution: "<< foot_pos << std::endl;



	Math::Vector3d	xe;
	Math::Matrix3d 	R;

	VectorNd q = VectorNd::Zero(model->q_size);

	forwardKinematics(q, xe, R, model, 5);

	std::cout << "xe : "<< xe << std::endl;
	std::cout << "R : "<< R << std::endl;
	Eigen::Vector3d xe_eig = xe;


	// NLOPT
	// Eigen::Vector3d target;
	// target << 0,0,0;
	// Data *data, d;
	// d.target = target;
	// d.body_id = 5;
	// d.model = model;

	// nlopt::opt opt(nlopt::LN_COBYLA, 6);

	// opt.set_min_objective(cost, NULL);
	// opt.set_xtol_rel(1e-4);
	// std::vector<double> q_nlopt(model->q_size);
	// std::fill(q_nlopt.begin(), q_nlopt.end(), 0);


	// double minf;
	// // Run the optimizer
	// t_start = std::chrono::high_resolution_clock::now();
	// std::cout << opt.optimize(q_nlopt, minf) << std::endl;
	// t_end = std::chrono::high_resolution_clock::now();
	// elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	// std::cout << "Time taken for IK nlopt [ms] : " << elapsed_time_ms  << std::endl;

	// // Print the results
	// for(int i = 0; i < q_nlopt.size(); i++) {
	// 	std::cout << "q indx : " << i << " : " << q_nlopt[i] << std::endl;
	// }

	delete model;
	return 0;
}





// int main() {
//   nlopt::opt opt(nlopt::LN_COBYLA, 6); //Providing gradients will decrease computational time required
//   opt.set_min_objective(cost, NULL);
//   opt.set_xtol_rel(1e-4);
//   std::vector<double> x(6);
//   x[0] = 0.0; x[1] = 0.0; x[2] = 0.0; x[3] = 0.0; x[4] = 0.0; x[5] = 0.0; x[6] = 0.0;
//   double minf;
//   // Run the optimizer
//   auto start = std::chrono::high_resolution_clock::now();
//   std::cout << opt.optimize(x, minf) << std::endl;
//   auto stop = std::chrono::high_resolution_clock::now();
//   auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
//   std::cout << "Time taken: " << duration.count() << std::endl;
//   // Print the results
//   std::cout << "Iteration count: " << count << std::endl;
//   std::cout << "Found minimum at f(" << x[0] << "," << x[1] << ","<< x[2] << "," << x[3] << "," << x[4] << "," << x[5] << "," << x[5] << ") = "
//             << std::setprecision(10) << minf <<std::endl;
// }

