#ifndef PINHOLE_H
#define PINHOLE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <Eigen/Geometry>

#include "../common/GeometricTransf.h"
#include "../common/StringUtils.h"

namespace at 
{
  class Pinhole
  {
 
  // TODO: Constructors using CAHV parameters need to extract the camera intrinsics and the extrinsic correctly
  // Rotation extraction is implemented but untested

  public:

    Pinhole(){};

    // Load the parameters from a file
    // This will throw exceptions if it fails
    // Object is deleted after throwing an exception
    // 0 means that the file did not open
    // 1 means that the file was misconfigured
    Pinhole(std::string cameraCalibrationFilename);  
    // All parameters, intrinsic and extrinsic, are specified
    Pinhole(int cols, int rows, double optical_center[2], double fx, double fy,
	    double rotation[9], double position[3]);
    // Use no rotation (identity matrix) and location at origin if no extrinsics are sent
    Pinhole(int cols, int rows, double optical_center[2], double fx, double fy);

   // Uses the parameters of a cahv model to construct a pinhole model
    Pinhole (double c[3], double a[3], double h[3], double v[3], int cols, int rows);
    Pinhole (double c[3], double a[3], double h[3], double v[3], int cols, int rows, double pixelSize);
  
    // Currently, nothing to destruct
    ~Pinhole(){};

    // Writes a file with the camera parameters that can be loaded later
    // Will overwrite any file of the same name
    void writeFile(std::string outputFilename);
    
    // Sets new parameters
    void setNew(int cols, int rows, double optical_center[2], double fx, double fy,
		double rotation[9], double position[3]);
 

    //Given a pixel, this function returns a ray that eminates from the camera center
    //and passes through that pixel
    //The length of the ray is determined by lambda, which will default to 100 if no lambda is sent
    //The camera position is also taken into account when calculating the ray.
    //The equation is as follows: end_point = position + lambda * simple_vector
    std::vector<float> pixelToVector(std::vector<int> pixel, double lambda);
    std::vector<float> pixelToVector(std::vector<int> pixel);

    //Given a position within the world, this function will return the pixel that corresponds to that
    //world location
    // WARNING: this function has not been tested extensively, it may have bugs
    std::vector<float> vectorToPixel(std::vector<float> v);

    void print();

    //get functions
    //intrinsic params
    double  GetFieldOfViewHor(){return m_fovX;};
    double  GetFieldOfViewVer(){return m_fovY;};
    double  GetFocalLengthHor(){return m_fx;};
    double  GetFocalLengthVer(){return m_fy;};
    double  GetOpticCenterHor(){return m_opticCenter[0];};
    double  GetOpticCenterVer(){return m_opticCenter[1];};
    bool    GetDistortion(std::vector<float>& get_distortion);
    int     GetRows () {return m_rows;};
    int     GetCols () {return m_cols;};
    void    GetIntrinsics(double get_intrinsics[3][3]);
    void    GetModel(double get_optical_center[2], double* get_fov_x, double* get_fov_y);

    //extrinsics
    double  GetCameraPositionX () {return m_position[0];};
    double  GetCameraPositionY () {return m_position[1];};
    double  GetCameraPositionZ () {return m_position[2];};
    void    GetRotation (double get_rotation[3][3]);  

    //these two functions will be removed - START
    void    GetSize(int& get_cols, int& get_rows);
    void    GetPosition(std::vector<float>& get_position);
    //these two functions will be removed - END
    
    // mosaic_processing needs this function
    // Given cahv parameters, this function will return the fov in radians
    void    GetPinholeModel(double C[3], double A[3], double H[3], double V[3], int cols, int rows,
			   float pixelSize, double optical_center[2], double *fov_X, double *fov_Y);

  private:

    int    m_cols, m_rows;
    float  m_pixelSize; 
    double m_fovX, m_fovY;
    double m_fx, m_fy;
    double m_opticCenter[2];
    std::vector<float> m_distortion;
    Eigen::Matrix3f m_camIntrinsics;
   
    std::vector<float> m_position;
    Eigen::Matrix3f m_rotation;
   

    // Send to common utility?
    double dotprod3(double v1[3], double v2[3]);
    void   crossprod3(double v1[3], double v2[3], double result[3]);
    double magnitude(double v[3]);
    void   scalar_mult(double s, double v1[3], double result[3]);
    void   add_vec3(double v1[3], double v2[3], double result[3]);
    void   sub_vec3(double v1[3], double v2[3], double result[3]);


    
  };
}
#endif
