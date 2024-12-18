#ifndef GEOMETRIC_TRANSF_H_
#define GEOMETRIC_TRANSF_H_

// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>

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

    void quaternionToRotation( double quaternion[4], double matrix[3][3] );

    // code from www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    // calculates a quaternion from a rotation matrix
    // Assumes the input matrix is special orthogonal
    void rotationToQuaternion( double a[3][3], double quaternion[4] );

    //rotates the input vector by the quaternion
    void rotateWithQuaternion( double input[3], double output[3], double quaternion[4] );

    Eigen::Matrix3f computeRotationFromXYZRPY(std::vector<float> xyzrpy);
    Eigen::Matrix4f xyzrpyToGeneralTransform(std::vector<float> xyzrpy);

    Eigen::Vector3d rotationMatrixToEulerXyz(Eigen::Matrix3d const& rotationMatrix);
    Eigen::Vector3f rotationMatrixToEulerXyzFloat(Eigen::Matrix3f const& rotationMatrix);
    void generalTransformToXYZRPY(std::vector<float> &xyzrpy, const Eigen::Matrix4f &xform);
    std::vector<float> transformXYZRPY(std::vector<float> xyzrpySrc, Eigen::Matrix4f transform);
}

#endif	// GEOMETRIC_TRANSF_H_
