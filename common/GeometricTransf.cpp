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
#include <Eigen/Dense>

static const double S_EPSILON = 0.0000001;

using namespace std;

namespace at
{
    // calculations from euclideanspace.com
    // function to take a quaternion and compute the 3x3 rotation matrix
    void quaternionToRotation( double quaternion[4], double matrix[3][3] )
    {
        //normalize
        double mag = sqrt(quaternion[0]*quaternion[0]+
			  quaternion[1]*quaternion[1]+
			  quaternion[2]*quaternion[2]+
			  quaternion[3]*quaternion[3]);

        double w = quaternion[0]/mag;
	double x = quaternion[1]/mag;
	double y = quaternion[2]/mag;
	double z = quaternion[3]/mag;

        double xx = x*x, xy = x*y, xz = x*z, xw = x*w;
        double yy = y*y, yz = y*z, yw = y*w;
        double zz = z*z, zw = z*w;

        matrix[0][0] = 1 - 2 * ( yy + zz );
        matrix[0][1] = 2 * ( xy - zw );
        matrix[0][2] = 2 * ( xz + yw );

        matrix[1][0] = 2 * ( xy + zw );
        matrix[1][1] = 1 - 2 * ( xx + zz );
        matrix[1][2] = 2 * ( yz - xw );

        matrix[2][0] = 2 * ( xz - yw );
        matrix[2][1] = 2 * ( yz + xw );
        matrix[2][2] = 1 - 2 * ( xx + yy );

    }

    void quaternionToMatrix(double quaternion[4], double m[3][3])
    {
        // quaternion for rotating the terrain into the local level frame
        double q0 = quaternion[0];	// w
        double q1 = quaternion[1];	// x
        double q2 = quaternion[2];	// y
        double q3 = quaternion[3];	// z
        static const double kQNormThreshold = 1.0e-3;
        double qNorm, q02, q0sqr2, q0q12, q0q22, q0q32, q1q22, q1q32, q2q32;

        // compute the unit vector k^ = Q / |Q|
        qNorm = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
        if (kQNormThreshold < qNorm)
        {
            q0 /= qNorm;
            q1 /= qNorm;
            q2 /= qNorm;
            q3 /= qNorm;

            // Paul Henning's Code to compute the equivalent Matrix
            // coeficients
            q02 = 2.0 * q0;
            q0sqr2 = q02 * q0;
            q0q12 = (q1 == 0.0) ? 0.0 : q02 * q1;
            q0q22 = (q2 == 0.0) ? 0.0 : q02 * q2;
            q0q32 = (q3 == 0.0) ? 0.0 : q02 * q3;
            q1q22 = 2.0*q1*q2;
            q1q32 = 2.0*q1*q3;
            q2q32 = 2.0*q2*q3;
    
            // This is for a matrix consistent with pre-multiplication by
            // row-vectors
            /*
              m[0][0] = q0sqr2 + q1*q1*2.0 - 1.0;	// 2*(q0^2 + q1^2) - 1
              m[0][1] = q1q22 + q0q32;	// = 2*q1*q2 + 2*q0*q3 = 2(x*y + w*z)
              m[0][2] = q1q32 - q0q22;	// = 2*q1*q3 - 2*q0*q2) = 2(x*z - w*y)

              m[1][0] = q1q22 - q0q32;	// 2.0*(q1*q2 - q0*q3) = 2*(x*y - w*z)
              m[1][1] = q0sqr2 + q2*q2*2.0 - 1.0;	// 2*(q0^2 + q2^2) - 1
              m[1][2] = q2q32 + q0q12;

              m[2][0] = q1q32 + q0q22;	// = 2*(q1*q3 + q0*q2) = 2*(x*z + w*y)
              m[2][1] = q2q32 - q0q12;	// = 2*(q2*q3 - q0*q1) = 2*(y*z - w*x)
              m[2][2] = q0sqr2 + q3*q3*2.0 - 1.0;	// 2*(q0^2 + q3^2) - 1
            */

            //this is its transpose
            m[0][0] = q0sqr2 + q1*q1*2.0 - 1.0;	// 2*(q0^2 + q1^2) - 1
            m[1][0] = q1q22 + q0q32;	// = 2*q1*q2 + 2*q0*q3 = 2(x*y + w*z)
            m[2][0] = q1q32 - q0q22;	// = 2*q1*q3 - 2*q0*q2) = 2(x*z - w*y)

            m[0][1] = q1q22 - q0q32;	// 2.0*(q1*q2 - q0*q3) = 2*(x*y - w*z)
            m[1][1] = q0sqr2 + q2*q2*2.0 - 1.0;	// 2*(q0^2 + q2^2) - 1
            m[2][1] = q2q32 + q0q12;

            m[0][2] = q1q32 + q0q22;	// = 2*(q1*q3 + q0*q2) = 2*(x*z + w*y)
            m[1][2] = q2q32 - q0q12;	// = 2*(q2*q3 - q0*q1) = 2*(y*z - w*x)
            m[2][2] = q0sqr2 + q3*q3*2.0 - 1.0;	// 2*(q0^2 + q3^2) - 1
	    /*
            cout<<"++++++++++++ROTATION MARIX++++++++++++++"<<endl;
            cout<<m[0][0]<<" "<<m[0][1]<<" "<<m[0][2]<<endl;
            cout<<m[1][0]<<" "<<m[1][1]<<" "<<m[1][2]<<endl;
            cout<<m[2][0]<<" "<<m[2][1]<<" "<<m[2][2]<<endl;
	    */
        }
        else				// the quaternion is not valid
        {
            printf("ERROR in QuaternionToMatrix: quaternion is not valid."
                   " Aborting.\n");
        }
    }

    // code from www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    // calculates a quaternion from a rotation matrix
    // Assumes the input matrix is special orthogonal
    void rotationToQuaternion( double a[3][3], double quaternion[4] )
    {

        double trace = a[0][0] + a[1][1] + a[2][2];
        if( trace > 0 ) {
            double s = 0.5 / sqrt(trace+ 1.0);
            quaternion[0] = 0.25 / s;
            quaternion[1] = ( a[2][1] - a[1][2] ) * s;
            quaternion[2] = ( a[0][2] - a[2][0] ) * s;
            quaternion[3] = ( a[1][0] - a[0][1] ) * s;
        } else {
            if ( a[0][0] > a[1][1] && a[0][0] > a[2][2] ) {
                double s = 2.0 * sqrt( 1.0 + a[0][0] - a[1][1] - a[2][2]);
                quaternion[0] = (a[2][1] - a[1][2] ) / s;
                quaternion[1] = 0.25 * s;
                quaternion[2] = (a[0][1] + a[1][0] ) / s;
                quaternion[3] = (a[0][2] + a[2][0] ) / s;
            } else if (a[1][1] > a[2][2]) {
                double s = 2.0 * sqrt( 1.0 + a[1][1] - a[0][0] - a[2][2]);
                quaternion[0] = (a[0][2] - a[2][0] ) / s;
                quaternion[1] = (a[0][1] + a[1][0] ) / s;
                quaternion[2] = 0.25 * s;
                quaternion[3] = (a[1][2] + a[2][1] ) / s;
            } else {
                double s = 2.0 * sqrt( 1.0 + a[2][2] - a[0][0] - a[1][1] );
                quaternion[0] = (a[1][0] - a[0][1] ) / s;
                quaternion[1] = (a[0][2] + a[2][0] ) / s;
                quaternion[2] = (a[1][2] + a[2][1] ) / s;
                quaternion[3] = 0.25 * s;
            }
        }

    }


    //rotates the input vector by the quaternion
    void rotateWithQuaternion( double input[3], double output[3], double quaternion[4] )
    {

        double M[3][3];
        //quaternionToRotation( quaternion, M ); 
        quaternionToMatrix(quaternion, M);
        double a0 = M[0][0]*input[0] + M[0][1]*input[1] + M[0][2]*input[2];
        double a1 = M[1][0]*input[0] + M[1][1]*input[1] + M[1][2]*input[2];
        double a2 = M[2][0]*input[0] + M[2][1]*input[1] + M[2][2]*input[2];

        output[0] = a0;
        output[1] = a1;
        output[2] = a2;

    }
    
  Eigen::Matrix3f computeRotationFromXYZRPY(vector<float> xyzrpy)
  {
    //(R_x*R_y*R_z)v
    //(r_roll*r_pitch*r_yaw)v
    //add the odometry offset to all points
    float roll =xyzrpy[3];
    float pitch=xyzrpy[4];
    float yaw  =xyzrpy[5];
    Eigen::Matrix3f r_roll; 
    Eigen::Matrix3f r_pitch;
    Eigen::Matrix3f r_yaw;
    Eigen::Matrix3f rot;
    r_roll(0,0)=1;
    r_roll(0,1)=0;
    r_roll(0,2)=0;
    r_roll(1,0)=0;
    r_roll(1,1)=cos(roll);
    r_roll(1,2)=-sin(roll);
    r_roll(2,0)=0;
    r_roll(2,1)=sin(roll);
    r_roll(2,2)=cos(roll);
    
    r_pitch(0,0)=cos(pitch);
    r_pitch(0,1)=0;
    r_pitch(0,2)=sin(pitch);
    r_pitch(1,0)=0;
    r_pitch(1,1)=1;
    r_pitch(1,2)=0;
    r_pitch(2,0)=-sin(pitch);
    r_pitch(2,1)=0;
    r_pitch(2,2)=cos(pitch);
    
    r_yaw(0,0)=cos(yaw);
    r_yaw(0,1)=-sin(yaw);
    r_yaw(0,2)=0;
    r_yaw(1,0)=sin(yaw);
    r_yaw(1,1)=cos(yaw);
    r_yaw(1,2)=0;
    r_yaw(2,0)=0;
    r_yaw(2,1)=0;
    r_yaw(2,2)=1;
    
//this is a valid convention but is the X-Y-Z *local* (aka Z-Y-X *global*) one
//that's not what roversw/vw uses nor what the rotation->euler function below
//uses
//    rot=r_roll*r_pitch*r_yaw;

    rot=r_yaw*r_pitch*r_roll;
    
    return rot;
  }
   
  /*
    Eigen::Matrix4f 
    HorizonMapGeo::rollPitchYawToRotationMatrix(
        const double& roll, const double& pitch, const double& yaw)
    {
        double sy = sin(yaw);
        double cy = cos(yaw);
	Eigen::Matrix4f yawMatrix;
	yawMatrix << 
	  cy,  -sy,   0,   0,
	  sy,   cy,   0,   0,
  	  0,    0,    1,   0,
	  0,    0,    0,   1;

        double sp = sin(pitch);
        double cp = cos(pitch);
	Eigen::Matrix4f pitchMatrix;
	pitchMatrix <<
	  cp,    0,    sp,  0,
	  0,     1,    0,   0,
	  -sp,   0,    cp,  0,
	  0,     0,    0,   1;

        double sr = sin(roll);
        double cr = cos(roll);
	Eigen::Matrix4f rollMatrix;
	rollMatrix <<
	  1,    0,     0,   0,
	  0,   cr,   -sr,   0,
	  0,   sr,    cr,   0,
          0,    0,     0,   1;

        return yawMatrix * pitchMatrix * rollMatrix;
    }

  */


  Eigen::Matrix4f xyzrpyToGeneralTransform(std::vector<float> xyzrpy)
  {
    Eigen::Matrix4f T;
    //Eigen::Matrix3f rotation = computeRotation(xyzrpy);
    Eigen::Matrix3f rotation = computeRotationFromXYZRPY(xyzrpy);
    T(0,0)=rotation(0,0); T(0,1)=rotation(0,1); T(0,2)=rotation(0,2); T(0,3)=xyzrpy[0]; 
    T(1,0)=rotation(1,0); T(1,1)=rotation(1,1); T(1,2)=rotation(1,2); T(1,3)=xyzrpy[1]; 
    T(2,0)=rotation(2,0); T(2,1)=rotation(2,1); T(2,2)=rotation(2,2); T(2,3)=xyzrpy[2];
    T(3,0)=0;             T(3,1)=0;             T(3,2)=0;             T(3,3)=1;
    return T;
  }

  Eigen::Vector3d rotationMatrixToEulerXyz(Eigen::Matrix3d const& rotationMatrix)
  {
    double phi = 0.;
    double const omega = asin(-rotationMatrix(2,0));
    double kappa;

    if (fabs(omega - M_PI_2) < S_EPSILON) {
      kappa = atan2( rotationMatrix(1,2),  rotationMatrix(0,2));
    }
    if (fabs(omega + M_PI_2) < S_EPSILON ) {
      kappa = atan2(-rotationMatrix(1,2), -rotationMatrix(0,2));
    }
    else {
      phi   = atan2(rotationMatrix(2,1), rotationMatrix(2,2));
      kappa = atan2(rotationMatrix(1,0), rotationMatrix(0,0));
    }
    return Eigen::Vector3d(phi, omega, kappa);
  }

  Eigen::Vector3f rotationMatrixToEulerXyzFloat(Eigen::Matrix3f const& rotationMatrix)
  {
    float phi = 0.;
    float const omega = asin(-rotationMatrix(2,0));
    float kappa;

    if (fabs(omega - M_PI_2) < S_EPSILON) {
      kappa = atan2( rotationMatrix(1,2),  rotationMatrix(0,2));
    }
    if (fabs(omega + M_PI_2) < S_EPSILON ) {
      kappa = atan2(-rotationMatrix(1,2), -rotationMatrix(0,2));
    }
    else {
      phi   = atan2(rotationMatrix(2,1), rotationMatrix(2,2));
      kappa = atan2(rotationMatrix(1,0), rotationMatrix(0,0));
    }
    return Eigen::Vector3f(phi, omega, kappa);
  }

  void generalTransformToXYZRPY(std::vector<float> &xyzrpy, const Eigen::Matrix4f &xform)
  {
    Eigen::Vector3f xyz = xform.block(0,3,3,1);
    Eigen::Vector3f rpy = rotationMatrixToEulerXyz(xform.block(0,0,3,3).cast<double>()).cast<float>();

    xyzrpy.resize(6);
    xyzrpy[0] = xyz.x();
    xyzrpy[1] = xyz.y();
    xyzrpy[2] = xyz.z();
    xyzrpy[3] = rpy.x();
    xyzrpy[4] = rpy.y();
    xyzrpy[5] = rpy.z();
  }

  vector<float> transformXYZRPY(vector<float> xyzrpySrc, Eigen::Matrix4f transform)
  {
    vector<float> xyzrpyDst;
    xyzrpyDst.resize(6);
    
    Eigen::Matrix4f transformSrc  = xyzrpyToGeneralTransform(xyzrpySrc);
    Eigen::Matrix4f transformDst = transform*transformSrc;
    generalTransformToXYZRPY(xyzrpyDst, transformDst);
    
    return xyzrpyDst;
  }
 

}
