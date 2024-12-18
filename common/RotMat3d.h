#ifndef COMMON_ROTMAT3D_H
#define COMMON_ROTMAT3D_H

#include <Eigen/Dense>

/**
 **** EULER ANGLE CONVENTION ****
 *
 * Unless otherwise specified, euler angles are in these equivalent conventions:
 *
 *  - X-Y-Z fixed
 *    Rotate around the original X axis, then the original Y axis, then the
 *    original Z axis.
 *
 *  - Z-Y-X local
 *    Rotate around the original Z axis, then the new Y axis, then the new-new
 *    X axis.
 *
 *  Angles are provided as (theta_x,theta_y,theta_z) = (roll,pitch,yaw) in rads
 */

namespace at
{

class RotMat3d
{
  public:
    //rotation matrices corresponding to a single rotation around x, y, or z:
    //R_{x,y,z}
    static Eigen::Matrix3d Rx(double theta);
    static Eigen::Matrix3d Ry(double theta);
    static Eigen::Matrix3d Rz(double theta);

    static void Rx(Eigen::Matrix3d &out, double theta);
    static void Ry(Eigen::Matrix3d &out, double theta);
    static void Rz(Eigen::Matrix3d &out, double theta);

    //partial derivative matrix corresponding to
    //dR_{x,y,z}/dtheta
    //for a single rotation around x, y, or z
    static Eigen::Matrix3d dRx_dtheta(double theta);
    static Eigen::Matrix3d dRy_dtheta(double theta);
    static Eigen::Matrix3d dRz_dtheta(double theta);

    static void dRx_dtheta(Eigen::Matrix3d &out, double theta);
    static void dRy_dtheta(Eigen::Matrix3d &out, double theta);
    static void dRz_dtheta(Eigen::Matrix3d &out, double theta);

    //partial derivative matrix corresponding to
    //dR(theta[x,y,z])/dtheta_{x,y, or z}
    //at an arbitrary current rotation specified by input euler vector theta
    //  dR_dthetax: R_z * R_y * dR_x/dtheta_x
    //  dR_dthetay: R_z * dR_y/dtheta_y * R_x
    //  dR_dthetaz: dR_z/dtheta_z * R_y * R_x
    static Eigen::Matrix3d dR_dthetax(const Eigen::Vector3d &theta);
    static Eigen::Matrix3d dR_dthetay(const Eigen::Vector3d &theta);
    static Eigen::Matrix3d dR_dthetaz(const Eigen::Vector3d &theta);

    static void dR_dthetax(Eigen::Matrix3d &out, const Eigen::Vector3d &theta);
    static void dR_dthetay(Eigen::Matrix3d &out, const Eigen::Vector3d &theta);
    static void dR_dthetaz(Eigen::Matrix3d &out, const Eigen::Vector3d &theta);

    //[rate of] change matrix corresponding to
    //dR(theta[x,y,z])/dtheta_x * thetadot_x
    // + dR(theta[x,y,z])/dtheta_y * thetadot_y
    // + dR(theta[x,y,z])/dtheta_z * thetadot_z
    //at an arbitrary current rotation specified by input euler vector theta
    //and a current angular rate specified by input vector thetadot
    static Eigen::Matrix3d dR_dtheta(const Eigen::Vector3d &theta, const Eigen::Vector3d &theta_dot);
    static void dR_dtheta(Eigen::Matrix3d &out, const Eigen::Vector3d &theta, const Eigen::Vector3d &theta_dot);

    //jacobian matrix of a rotation represented as a rotation matrix R
    //  whose output is represented in euler angles
    //  with respect to an euler angle change in R
    //  given the supplied derivative matrices of R wrt x, y, and z euler angles
    //corresponding to
    //d(rotation2euler(R))/dtheta_[x,y,z] | dR/dtheta_[x,y,z]
    static Eigen::Matrix3d dthetaR_dtheta(const Eigen::Matrix3d &R, const Eigen::Matrix3d &dR_dthetax, const Eigen::Matrix3d &dR_dthetay, const Eigen::Matrix3d &dR_dthetaz);
    static void dthetaR_dtheta(Eigen::Matrix3d &out, const Eigen::Matrix3d &R, const Eigen::Matrix3d &dR_dthetax, const Eigen::Matrix3d &dR_dthetay, const Eigen::Matrix3d &dR_dthetaz);

};

}; /* namespace at */

#endif //COMMON_ROTMAT3D_H

