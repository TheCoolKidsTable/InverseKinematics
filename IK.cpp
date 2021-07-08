#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <nlopt.hpp>
#include <vector>
#include <numeric>
#include <Eigen/Geometry>
#include <chrono>

int count;

Eigen::Matrix3d skew(Eigen::Vector3d u) {
  Eigen::Matrix3d skew;
  skew << 0, -u(2), u(1),
      u(2), 0, -u(0),
      -u(1), u(0), 0;
  return skew;
}

// Generator matrix
Eigen::Matrix4d generatorMatrix(Eigen::MatrixXd x) {
  Eigen::Vector3d x13;
  x13 << x(1), x(2), x(3);
  Eigen::Vector3d x46;
  x46 << x(4), x(5), x(6);
  Eigen::Matrix3d skew_x46 = skew(x46);
  Eigen::Matrix4d GG;
  GG << skew_x46(1), skew_x46(2), skew_x46(3), x13(1),
       skew_x46(4), skew_x46(5), skew_x46(6), x13(2),
       skew_x46(7), skew_x46(8), skew_x46(9), x13(3),
       0.0, 0.0, 0.0, 0.0;
  return GG;
}

double cost(const std::vector<double> &q, std::vector<double> &grad, void *data) {
        Eigen::Affine3d runningTransform = Eigen::Affine3d::Identity();
        double hip_yaw = 0.055;
        double hip_roll = 0.06;
        double hip_pitch = 0;
        double upper_leg = 0.2;
        double lower_leg = 0.2;
        double ankle = 0;
        int negativeIfRight = -1;
        Eigen::MatrixXd twist_vector(6,1);
        twist_vector << 0,0,0,0,0,1;
        Eigen::Matrix4d derivativeMatrix;
        derivativeMatrix = generatorMatrix(twist_vector);
        // Hip yaw
        Eigen::Affine3d H01 = Eigen::Affine3d::Identity();
        H01 = H01.translate(Eigen::Vector3d(0.0, negativeIfRight*hip_yaw, 0.0));
        // Calculate derivative with respect to q1
        Eigen::Matrix4d dH01dq1 = Eigen::Matrix4d::Zero();
        // Hip roll
        Eigen::Affine3d H12 = Eigen::Affine3d::Identity();
        H12 = H12.rotate(Eigen::AngleAxisd(q[0], Eigen::Vector3d::UnitZ()));
        H12 = H12.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitY()));
        H12 = H12.translate(Eigen::Vector3d(hip_roll, 0.0, 0.0));
        // Calculate derivative with respect to q2
        Eigen::Matrix4d dH12dq2;
        dH12dq2 = derivativeMatrix*H12.matrix();
        // Hip pitch
        Eigen::Affine3d H23 = Eigen::Affine3d::Identity();
        H23 = H23.rotate(Eigen::AngleAxisd(q[1], Eigen::Vector3d::UnitZ()));
        H23 = H23.rotate(Eigen::AngleAxisd(-M_PI_2, Eigen::Vector3d::UnitX()));
        H23 = H23.translate(Eigen::Vector3d(hip_pitch, 0.0, 0.0));
        // Calculate derivative with respect to q3
        Eigen::Matrix4d dH23dq3;
        dH23dq3 = derivativeMatrix*H23.matrix();
        // Upper leg
        Eigen::Affine3d H34 = Eigen::Affine3d::Identity();
        H34 = H34.rotate(Eigen::AngleAxisd(q[2], Eigen::Vector3d::UnitZ()));
        H34 = H34.translate(Eigen::Vector3d(upper_leg, 0.0, 0.0));
        // Calculate derivative with respect to q4
        Eigen::Matrix4d dH34dq4;
        dH34dq4 = derivativeMatrix*H34.matrix();
        // Lower leg
        Eigen::Affine3d H45 = Eigen::Affine3d::Identity();
        H45 = H45.rotate(Eigen::AngleAxisd(q[3], Eigen::Vector3d::UnitZ()));
        H45 = H45.translate(Eigen::Vector3d(lower_leg, 0.0, 0.0));
        // Calculate derivative with respect to q5
        Eigen::Matrix4d dH45dq5;
        dH45dq5 = derivativeMatrix*H45.matrix();
        // Ankle
        Eigen::Affine3d H56 = Eigen::Affine3d::Identity();
        H56 = H56.rotate(Eigen::AngleAxisd(q[4], Eigen::Vector3d::UnitZ()));
        H56 = H56.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitX()));
        H56 = H56.translate(Eigen::Vector3d(ankle, 0.0, 0.0));
        // Calculate derivative with respect to q6
        Eigen::Matrix4d dH56dq6;
        dH56dq6 = derivativeMatrix*H56.matrix();
        // Foot
        Eigen::Affine3d H6R = Eigen::Affine3d::Identity();
        H6R = H6R.rotate(Eigen::AngleAxisd(q[5], Eigen::Vector3d::UnitZ()));
        // Calculate derivative with respect to q7
        Eigen::Matrix4d dH6Rdq7;
        dH6Rdq7 = derivativeMatrix*H6R.matrix();

        // Use product rule to calculate analytical jacobian
        Eigen::Matrix4d dH0Rdq1 = dH01dq1*H12.matrix()*H23.matrix()*H34.matrix()*H45.matrix()*H56.matrix()*H6R.matrix();
        Eigen::Matrix4d dH0Rdq2 = H01.matrix()*dH12dq2*H23.matrix()*H34.matrix()*H45.matrix()*H56.matrix()*H6R.matrix();
        Eigen::Matrix4d dH0Rdq3 = H01.matrix()*H12.matrix()*dH23dq3*H34.matrix()*H45.matrix()*H56.matrix()*H6R.matrix();
        Eigen::Matrix4d dH0Rdq4 = H01.matrix()*H12.matrix()*H23.matrix()*dH34dq4*H45.matrix()*H56.matrix()*H6R.matrix();
        Eigen::Matrix4d dH0Rdq5 = H01.matrix()*H12.matrix()*H23.matrix()*H34.matrix()*dH45dq5*H56.matrix()*H6R.matrix();
        Eigen::Matrix4d dH0Rdq6 = H01.matrix()*H12.matrix()*H23.matrix()*H34.matrix()*H45.matrix()*dH56dq6*H6R.matrix();
        Eigen::Matrix4d dH0Rdq7 = H01.matrix()*H12.matrix()*H23.matrix()*H34.matrix()*H45.matrix()*H56.matrix()*dH6Rdq7.matrix();
        Eigen::MatrixXd JA(56,4);
        JA << dH0Rdq1,dH0Rdq2,dH0Rdq3,dH0Rdq4,dH0Rdq5,dH0Rdq6,dH0Rdq7;
        Eigen::Affine3d H0R = Eigen::Affine3d::Identity();
        H0R = H01*H12*H23*H34*H45*H56*H6R;

        count++;
        Eigen::Vector3d xe = H0R.translation();
        Eigen::Vector3d xd(0.4,-0.055,0);
        Eigen::Vector3d cost = xe - xd;
        double K = 1;
        // Gradient
        if(!grad.empty()){
        grad[0] = 2*dH0Rdq1.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq1.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq1.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[1] = 2*dH0Rdq2.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq2.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq2.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[2] = 2*dH0Rdq3.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq3.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq3.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[3] = 2*dH0Rdq4.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq4.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq4.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[4] = 2*dH0Rdq5.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq5.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq5.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[5] = 2*dH0Rdq6.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq6.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq6.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        grad[6] = 2*dH0Rdq7.coeff(1,4)*(xe.coeff(1,1)-xd.coeff(1,1))+2*dH0Rdq7.coeff(2,4)*(xe.coeff(2,1)-xd.coeff(2,1))+2*dH0Rdq1.coeff(3,4)*(xe.coeff(3,1)-xd.coeff(3,1));
        }
        return cost.transpose()*K*cost;
}

int main() {
  nlopt::opt opt(nlopt::AUGLAG, 6); //Providing gradients will decrease computational time required
  opt.set_min_objective(cost, NULL);
  nlopt::opt local_opt(nlopt::LN_COBYLA, 6);
  opt.set_local_optimizer(local_opt);  
  opt.set_xtol_rel(1e-1);
  std::vector<double> x(6);
  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0; x[3] = 0.0; x[4] = 0.0; x[5] = 0.0; x[6] = 0.0;
  double minf;
  // Run the optimizer
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << opt.optimize(x, minf) << std::endl;
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken: " << duration.count() << std::endl;
  // Print the results
  std::cout << "Iteration count: " << count << std::endl;
  std::cout << "Found minimum at f(" << x[0] << "," << x[1] << ","<< x[2] << "," << x[3] << "," << x[4] << "," << x[5] << "," << x[5] << ") = "
            << std::setprecision(10) << minf <<std::endl;
}
