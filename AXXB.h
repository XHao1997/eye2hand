//
// Created by Hao on 25/06/2022.
//

#ifndef UNTITLED2_AXXB_H
#define UNTITLED2_AXXB_H

#endif //UNTITLED2_AXXB_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/SVD>

Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd &origin, const float er);

Eigen::Matrix3d Eye2HandCli(Eigen::MatrixXd A, Eigen::MatrixXd B);

Eigen::Quaterniond rotationMatrix2Quaterniond(Eigen::Matrix3d R);

Eigen::Matrix3d skew(Eigen::Quaternion<double> R);

Eigen::Matrix4d CalAA(Eigen::Quaterniond a, Eigen::Quaterniond b);


