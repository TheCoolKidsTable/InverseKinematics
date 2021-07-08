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



Eigen::Affine3d forwardKinematics(const std::vector<double> &q) {
        Eigen::Affine3d runningTransform = Eigen::Affine3d::Identity();
        double hip_yaw = 0.055;
        double hip_roll = 0.06;
        double hip_pitch = 0;
        double upper_leg = 0.2;
        double lower_leg = 0.2;
        double ankle = 0;
        int negativeIfRight = -1;
        // Hip yaw
        runningTransform = runningTransform.translate(Eigen::Vector3d(0.0, negativeIfRight*hip_yaw, 0.0));
        // Hip roll
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[0], Eigen::Vector3d::UnitZ()));
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitY()));
        runningTransform = runningTransform.translate(Eigen::Vector3d(hip_roll, 0.0, 0.0));
        // Hip pitch
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[1], Eigen::Vector3d::UnitZ()));
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(-M_PI_2, Eigen::Vector3d::UnitX()));
        runningTransform = runningTransform.translate(Eigen::Vector3d(hip_pitch, 0.0, 0.0));
        // Upper leg
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[2], Eigen::Vector3d::UnitZ()));
        runningTransform = runningTransform.translate(Eigen::Vector3d(upper_leg, 0.0, 0.0));
        // Lower leg
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[3], Eigen::Vector3d::UnitZ()));
        runningTransform = runningTransform.translate(Eigen::Vector3d(lower_leg, 0.0, 0.0));
        // Ankle
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[4], Eigen::Vector3d::UnitZ()));
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(M_PI_2, Eigen::Vector3d::UnitX()));
        runningTransform = runningTransform.translate(Eigen::Vector3d(ankle, 0.0, 0.0));
        // Foot
        runningTransform = runningTransform.rotate(Eigen::AngleAxisd(q[5], Eigen::Vector3d::UnitZ()));
        return runningTransform;
}

double cost(const std::vector<double> &q, std::vector<double> &grad, void *data)
{
  count++;
  Eigen::Vector3d xe = forwardKinematics(q).translation();
  Eigen::Vector3d xd(0.4,-0.055,0);
  Eigen::Vector3d cost = xe - xd;
  return cost.norm();
}

int main() {
  nlopt::opt opt(nlopt::LN_COBYLA, 6); //Providing gradients will decrease computational time required
  opt.set_min_objective(cost, NULL);
  opt.set_xtol_rel(1e-4);
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
