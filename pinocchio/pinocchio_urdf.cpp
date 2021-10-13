#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/center-of-mass.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/geometry.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/urdf/utils.hpp"
#include "pinocchio/parsers/urdf/model.hxx"

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
#include <math.h>
#include <nlopt.hpp>

#include <iostream>

typedef struct {
	pinocchio::Model model;
	Eigen::Affine3d target;
	int JOINT_ID;
  	pinocchio::Data data;
} params;



double cost(const std::vector<double> & q_nlopt, std::vector<double> &grad, void* data)
{
	params* p = reinterpret_cast<params*>(data);
	pinocchio::Model model = p->model;
	int JOINT_ID = p->JOINT_ID;
	pinocchio::Data pinocchio_data = p->data;

	Eigen::VectorXd q;
	q.resize(q_nlopt.size());
	for(int i = 0; i < q_nlopt.size(); i++) {
		q(i) = q_nlopt[i];
	}

	// calculate current end effector position
	pinocchio::forwardKinematics(model,pinocchio_data,q);
	Eigen::Vector3d xe = pinocchio_data.oMi[JOINT_ID].translation();
	Eigen::Matrix3d Re = pinocchio_data.oMi[JOINT_ID].rotation();
	std::cout << "R : " << Re << std::endl;
	Eigen::Quaternionf Qe(Re);

	//  x desired
	Eigen::Vector3d xd = p->target.translation();
	Eigen::MatrixXd position_cost = xe - xd;


	Eigen::MatrixXd orientation_cost = xe - xd;
	return position_cost.norm();
}



int main(int argc, char ** argv)
{
    using namespace pinocchio;

    // You should change here to set up your own URDF file or just pass it as an argument of this example.
    const std::string urdf_filename = argv[1];

    // Load the urdf model
    pinocchio::Model model;
    pinocchio::urdf::buildModel(urdf_filename,model, false);
    std::cout << "model name: " << model.name << std::endl;

    // Create data required by the algorithms
    pinocchio::Data data(model);

    // Sample a random configuration
    Eigen::VectorXd q = pinocchio::neutral(model);
	  q.setZero();
    std::cout << "q: " << q.transpose() << std::endl;
    // Perform the forward kinematics over the kinematic tree
    forwardKinematics(model,data,q);
    // Calculate Centre of Mass
    Eigen::MatrixXd CoM = pinocchio::centerOfMass(model,data,q,false);
    std::cout << "CoM : " << CoM << std::endl;
    // Print out the placement of each joint of the kinematic tree
    for(JointIndex joint_id = 0; joint_id < (JointIndex)model.njoints; ++joint_id)
      std::cout << std::setw(24) << std::left
                << model.names[joint_id] << ": "
                << std::fixed << std::setprecision(4)
                << data.oMi[joint_id].translation().transpose()
                << std::endl;

    // Inverse Kinematics
    const int JOINT_ID = 5;
    const pinocchio::SE3 oMdes(Eigen::Matrix3d::Identity(), Eigen::Vector3d(0.0294,  0.0533,  -0.4529));

    const double eps  = 1e-4;
    const int IT_MAX  = 1000;
    const double DT   = 1e-1;
    const double damp = 1e-6;

    pinocchio::Data::Matrix6x J(6,model.nv);
    J.setZero();

    bool success = false;
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    Vector6d err;
    Eigen::VectorXd v(model.nv);
    auto t_start_3 = std::chrono::high_resolution_clock::now();
    for (int i=0;;i++)
    {
      pinocchio::forwardKinematics(model,data,q);
      const pinocchio::SE3 dMi = oMdes.actInv(data.oMi[JOINT_ID]);
      err = pinocchio::log6(dMi).toVector();
      if(err.norm() < eps)
      {
        success = true;
        break;
      }
      if (i >= IT_MAX)
      {
        success = false;
        break;
      }
      pinocchio::computeJointJacobian(model,data,q,JOINT_ID,J);
      pinocchio::Data::Matrix6 JJt;
      JJt.noalias() = J * J.transpose();
      JJt.diagonal().array() += damp;
      v.noalias() = - J.transpose() * JJt.ldlt().solve(err);
      q = pinocchio::integrate(model,q,v*DT);
      if(!(i%10))
        std::cout << i << ": error = " << err.transpose() << std::endl;
    }

    auto t_end_3 = std::chrono::high_resolution_clock::now();
    auto elapsed_time_ms_3 = std::chrono::duration<double, std::milli>(t_end_3-t_start_3).count()/1000;
    std::cout << "Time taken for IK pinocchio [ms] : " << elapsed_time_ms_3  << std::endl;

    if(success)
      std::cout << "Convergence achieved!" << std::endl;
    else
      std::cout << "\nWarning: the iterative algorithm has not reached convergence to the desired precision" << std::endl;

    std::cout << "\nresult: " << q.transpose() << std::endl;
    std::cout << "\nfinal error: " << err.transpose() << std::endl;


  std::cout << "************************ NLopt Forward Kinematics ************************" << std::endl;

	Eigen::Affine3d target;
	Eigen::Vector3d orientation_target;
	Eigen::Vector3d position_target;
	position_target <<  0.0294,  0.0533,  -0.4;
	orientation_target << 0,1.571,0;
	target.translation() = position_target;
	target.linear() = Eigen::MatrixXd::Identity(3,3);


	params pinocchio_data[4] = {model, target, 5, data};

	nlopt::opt opt(nlopt::LN_COBYLA, 20);
	opt.set_min_objective(cost, &pinocchio_data[0]);
	opt.set_xtol_rel(1e-5);
	std::vector<double> q_nlopt(20);
  std::cout << "here" << std::endl;
	// Initial conditions
	std::fill(q_nlopt.begin(), q_nlopt.end(), 0);
	double minf;
	// Run the optimizer
	auto t_start_4 = std::chrono::high_resolution_clock::now();
	std::cout << opt.optimize(q_nlopt, minf) << std::endl;
	auto t_end_4 = std::chrono::high_resolution_clock::now();
	auto elapsed_time_ms_4 = std::chrono::duration<double, std::milli>(t_end_4-t_start_4).count();
	std::cout << "Time taken for IK nlopt [ms] : " << elapsed_time_ms_3  << std::endl;

	// Print the results
	for(int i = 0; i < q_nlopt.size(); i++) {
		std::cout << "q indx : " << i << " : " << q_nlopt[i] << std::endl;
	}

	q.resize(q_nlopt.size());
	for(int i = 0; i < q_nlopt.size(); i++) {
		q(i) = q_nlopt[i];
	}

  std::cout << "Result NLopt : " << std::endl;
  forwardKinematics(model,data,q);
  // Print out the placement of each joint of the kinematic tree
  for(JointIndex joint_id = 0; joint_id < (JointIndex)model.njoints; ++joint_id)
    std::cout << std::setw(24) << std::left
              << model.names[joint_id] << ": "
              << std::fixed << std::setprecision(4)
              << data.oMi[joint_id].translation().transpose()
              << std::endl;

}