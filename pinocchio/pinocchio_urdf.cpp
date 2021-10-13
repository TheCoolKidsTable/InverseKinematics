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

#include "rotation.hpp"
#include "matlab.hpp"

#include "mat.h"


typedef struct {
	pinocchio::Model model;
	Eigen::Affine3d target;
	int JOINT_ID;
  	pinocchio::Data data;
	std::vector<double> q0;
} params;



Eigen::Matrix3d skew (Eigen::Vector3d u) {
	Eigen::Matrix3d u_hat;
	u_hat << 0, -u(2), u(1),
		u(2), 0, -u(0),
		-u(1), u(0), 0;
		return u_hat;
}



double cost(const std::vector<double> & q_nlopt, std::vector<double> &grad, void* data)
{
	params* p = reinterpret_cast<params*>(data);
	pinocchio::Model model = p->model;
	int JOINT_ID = p->JOINT_ID;
	pinocchio::Data pinocchio_data = p->data;

	Eigen::VectorXd q;
	Eigen::VectorXd q0;
	q.resize(q_nlopt.size());
	q0.resize(q_nlopt.size());
	for(int i = 0; i < q_nlopt.size(); i++) {
		q(i) = q_nlopt[i];
		q0(i) = p->q0[i];
	}

	// calculate current end effector position
	pinocchio::forwardKinematics(model,pinocchio_data,q);

	// position cost
	Eigen::Vector3d xe = pinocchio_data.oMi[JOINT_ID].translation();
	Eigen::Vector3d xd = p->target.translation();
	Eigen::MatrixXd position_cost = xe - xd;

	// orientation cost
	Eigen::Matrix3d Re = pinocchio_data.oMi[JOINT_ID].rotation();
	Eigen::Quaterniond Qe(Re);
	Eigen::Quaterniond Qd(p->target.linear());

	// Robotics - Modelling, Planning and Control pg. 140
	// cost = ne*ed -ne*ee - S(ed)*ee
	double eta_e = Qe.w();
	double eta_d = Qd.w();
	Eigen::Vector3d epsilon_e(3,1);
	epsilon_e << Qe.x(), Qe.y(), Qe.z();
	Eigen::Vector3d epsilon_d(3,1);
	epsilon_d << Qd.x(), Qd.y(), Qd.z();
	Eigen::MatrixXd skew_of_epsilon_d(3,3);
	skew_of_epsilon_d = skew(epsilon_d);
	Eigen::Vector3d orientation_cost(3,1);
	orientation_cost = eta_e*epsilon_d - eta_d*epsilon_e - skew_of_epsilon_d*epsilon_e;

	// Initial condition cost
	Eigen::MatrixXd W(q_nlopt.size(),q_nlopt.size());
	W.diagonal().fill(0.0001);
	double initial_cost = (q - q0).transpose()*W*(q-q0);

	Eigen::MatrixXd K;
	K.diagonal().fill(0.01);


	return position_cost.norm() + orientation_cost.transpose()*orientation_cost + initial_cost;
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
		std::cout << std::setw(40) << std::left
				<< model.names[joint_id] 
				<< "JOINT_ID : " << joint_id
				<< " position : "
				<< std::fixed << std::setprecision(4)
				<< data.oMi[joint_id].translation().transpose()
				<< " orientation : " << data.oMi[joint_id].rotation().eulerAngles(2,1,0).transpose()
				<< std::endl;

    // Inverse Kinematics
    const int JOINT_ID = 6;
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

	// **** Position Target ****
	Eigen::Vector3d position_target;
	position_target <<  0.0294,  0.0523,  -0.3529;
	target.translation() = position_target;

	// **** Orientation Target ****
	// Target rotation matrix
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Rd;
	// Target euler angles
	Eigen::Matrix<double, Eigen::Dynamic, 1> Thetad;
	Thetad.resize(3,1);
	Thetad << 0, 1.571, 0.0000;
	rpy2rot(Thetad, Rd);
	target.linear() = Rd;

	// Initial conditions
	std::vector<double> q_nlopt(20);
	std::fill(q_nlopt.begin(), q_nlopt.end(), 0);
	params pinocchio_data[4] = {model, target, 6, data, q_nlopt};

	auto t_start_4 = std::chrono::high_resolution_clock::now();
	nlopt::opt opt(nlopt::LN_COBYLA, 20);
	opt.set_min_objective(cost, &pinocchio_data[0]);
	opt.set_xtol_rel(1e-4);
	double minf;
	// Run the optimizer
	std::cout << opt.optimize(q_nlopt, minf) << std::endl;
	auto t_end_4 = std::chrono::high_resolution_clock::now();
	auto elapsed_time_ms_4 = std::chrono::duration<double, std::milli>(t_end_4-t_start_4).count();
	std::cout << "Time taken for IK nlopt [ms] : " << elapsed_time_ms_4  << std::endl;

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
    std::cout << std::setw(40) << std::left
              << model.names[joint_id] 
			  << "JOINT_ID : " << joint_id
			  << " position : "
              << std::fixed << std::setprecision(4)
              << data.oMi[joint_id].translation().transpose()
			  << " orientation : " << data.oMi[joint_id].rotation().eulerAngles(2,1,0).transpose()
              << std::endl;

	std::cout << "Joint 6 position : " << data.oMi[6].translation().transpose() << std::endl;
	std::cout << "Joint 6 rotation matrix" << data.oMi[6].rotation() << std::endl;

	const std::string file_name = "../q.csv" ;
	const std::string vec_name = "vector of joint angles";
	std::ofstream file(file_name);
	save_vector_as_matrix(vec_name,q_nlopt,file);

	// system("visualize.sh");
	// std::cin.ignore();

}