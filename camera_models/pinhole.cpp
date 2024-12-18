#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include "pinhole.h"

using namespace std;

namespace at
{

    double
    Pinhole::dotprod3(double v1[3], double v2[3])
    {
        return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    }

    void
    Pinhole::crossprod3(double v1[3], double v2[3], double result[3])
    {
        // a x b = (a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1)
        result[0] = (v1[1] * v2[2] - v1[2] * v2[1]);
        result[1] = (v1[2] * v2[0] - v1[0] * v2[2]);
        result[2] = (v1[0] * v2[1] - v1[1] * v2[0]);
    }

    double
    Pinhole::magnitude(double v[3])
    {
        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    void
    Pinhole::scalar_mult(double s, double v1[3], double result[3])
    {
        result[0] = s * v1[0];
        result[1] = s * v1[1];
        result[2] = s * v1[2];
    }

    void
    Pinhole::add_vec3(double v1[3], double v2[3], double result[3])
    {
        result[0] = v1[0] + v2[0];
        result[1] = v1[1] + v2[1];
        result[2] = v1[2] + v2[2];
    }

    void 
    Pinhole::sub_vec3(double v1[3], double v2[3], double result[3])
    {
        result[0] = v1[0] - v2[0];
        result[1] = v1[1] - v2[1];
        result[2] = v1[2] - v2[2];
    }

 
    //constructor from file
    Pinhole::Pinhole(string cameraCalibrationFilename)
    {
      ifstream calibrationFile (cameraCalibrationFilename.c_str());
      string line;
      int found = 0;
      int require_parameters = 2;
      
      //set defalut parameters
      m_position.resize(3);
      
      m_position[0] = 0.0;
      m_position[1] = 0.0;
      m_position[2] = 0.0;
      
      m_rotation(0, 0) = 1.0;
      m_rotation(0, 1) = 0.0;
      m_rotation(0, 2) = 0.0;
      m_rotation(1, 0) = 0.0;
      m_rotation(1, 1) = 1.0;
      m_rotation(1, 2) = 0.0;
      m_rotation(2, 0) = 0.0;
      m_rotation(2, 1) = 0.0;
      m_rotation(2, 2) = 1.0;
      
      m_distortion.resize(4);
      m_distortion[0]=0.0;
      m_distortion[1]=0.0;
      m_distortion[2]=0.0;
      m_distortion[3]=0.0;
      
      if (calibrationFile.is_open()){
	while (!calibrationFile.eof()){
	  line.clear();
	  std::getline(calibrationFile, line);
	  
	  if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines
	    std::stringstream lineStream(line);
	    std::string keyword("");
	    lineStream >> keyword;
	    
	    if (keyword.compare(std::string("CAMERA_MATRIX"))==0){
	      found++;
	      
	      // Check number of parameters
	      vector<string> contents = split(line);
	      
	      if(contents.size() < 10){
		// Not enough parameters
		// Throw failed load
		throw 1;
	      }
	      
	      // Get the parameters
	      m_opticCenter[0] = atof(contents[3].c_str());
	      m_opticCenter[1] = atof(contents[6].c_str());
	      m_fx = atof(contents[1].c_str());
	      m_fy = atof(contents[5].c_str());
	    }
	    else if (keyword.compare(std::string("TRANSLATION"))==0){
	      
	      //found++;
	      
	      // Check number of parameters
	      vector<string> contents = split(line);
	      if(contents.size() < 4){
		// Not enough parameters
		// Throw failed load
		throw 1;
	      }
	      
	      // Get the parameters
	      m_position[0] = atof(contents[1].c_str());
	      m_position[1] = atof(contents[2].c_str());
	      m_position[2] = atof(contents[3].c_str());
	    }			
	    else if (keyword.compare(std::string("WIDTH_HEIGHT"))==0){
	      found++;
	      
	      // Check number of parameters
	      vector<string> contents = split(line);
	      if(contents.size() < 3){
		// Not enough parameters
		// Throw failed load
		throw 1;
	      }
	      
	      // Get the parameters
	      m_cols = atof(contents[1].c_str());
	      m_rows = atof(contents[2].c_str());
	    }
	    else if (keyword.compare(std::string("ROTATION"))==0){
	      //found++;
	      
	      // Check number of parameters
	      vector<string> contents = split(line);
	      
	      if(contents.size() < 10){
		// Not enough parameters
		// Throw failed load
		throw 1;
	      }
	      
	      // Get the parameters
	      m_rotation(0, 0) = atof(contents[1].c_str());
	      m_rotation(0, 1) = atof(contents[2].c_str());
	      m_rotation(0, 2) = atof(contents[3].c_str());
	      m_rotation(1, 0) = atof(contents[4].c_str());
	      m_rotation(1, 1) = atof(contents[5].c_str());
	      m_rotation(1, 2) = atof(contents[6].c_str());
	      m_rotation(2, 0) = atof(contents[7].c_str());
	      m_rotation(2, 1) = atof(contents[8].c_str());
	      m_rotation(2, 2) = atof(contents[9].c_str());
	    }
	    else if (keyword.compare(std::string("DISTORTION_COEFFICIENTS"))==0){
	      // WARNING: the distortion coefficients may have varying numbers, no check
	      //is being performed for number of coefficients
	      // Get the contents of the line
	      vector<string> contents = split(line);
	      m_distortion.resize(0);
	      // Copy all available coefficients
	      // First member of contents is the identifier, so ignore it
	      for(unsigned int i = 1; i < contents.size(); i++){
		m_distortion.push_back(atof(contents[i].c_str()));
	      }
	    }
	  }
	}
	calibrationFile.close();
	
	// Check for all parameters found
	if(found < require_parameters/* || (!found_distortion && use_distortion) */){
	  // Not enough parameters or distortion was not found when required
	  // Throw failed load
	  cout<<"Pinhole:: Required fields not found in pinhole file "<<cameraCalibrationFilename<<"."<<endl;
	  throw 1;
	}
	
	// Set the intrinsics
	m_camIntrinsics << m_fx, 0, m_opticCenter[0],
	  0, m_fy, m_opticCenter[1],
	  0, 0, 1;
      }
      else{
	cout<<"Pinhole:: Pinhole file "<<cameraCalibrationFilename<<"not found. Exiting."<<endl;
	// Throw failed open
	throw 0;
      }
    }
  
    Pinhole::Pinhole(int cols, int rows, double opticCenter[2], double fx, double fy)
    {
      m_position.resize(3);

      // Set the cols and rows
      m_cols = cols;
      m_rows = rows;

      // Set the optic center
      m_opticCenter[0]=opticCenter[0];
      m_opticCenter[1]=opticCenter[1];

      // Set the focal lengths
      m_fx = fx;
      m_fy = fy;

      // Set the rotation to identity
      m_rotation << 1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0;

      // Set the position
      m_position[0] = 0.0;
      m_position[1] = 0.0;
      m_position[2] = 0.0;

      // Set the intrinsics
      m_camIntrinsics << m_fx, 0,    m_opticCenter[0],
                         0,    m_fy, m_opticCenter[1],
                         0,    0,    1;
    };
  
  Pinhole::Pinhole(int cols, int rows, double opticCenter[2], double fx, double fy,
		   double rotation[9], double position[3])
  {
    
    m_position.resize(3);
    
    // Set the cols and rows
    m_cols = cols;
    m_rows = rows;
    
    // Set the optic center
    m_opticCenter[0]=opticCenter[0];
    m_opticCenter[1]=opticCenter[1];
    
    // Set the focal lengths
    m_fx = fx;
    m_fy = fy;
    
    // Set the rotation
    m_rotation << rotation[0], rotation[1], rotation[2],
                  rotation[3], rotation[4], rotation[5],
                  rotation[6], rotation[7], rotation[8];
    
    // Set the position
    m_position[0] = position[0];
    m_position[1] = position[1];
    m_position[2] = position[2];
    
    // Set the intrinsics
    m_camIntrinsics << m_fx, 0,    m_opticCenter[0],
                       0,    m_fy, m_opticCenter[1],
                       0,    0,    1;
  };
  
  //constructor from CAHV parameters
  Pinhole::Pinhole(double c[3], double a[3], double h[3], double v[3], int cols, int rows)
  {
    double x_axis[3] = {0.0, 0.0, 0.0}, y_axis[3] = {0.0, 0.0, 0.0};
    double AxH[3], AxV[3];
    crossprod3(a, h, AxH);
    crossprod3(a, v, AxV);
    double horizontal_scale = magnitude(AxH);
    double vertical_scale = magnitude(AxV);
    m_opticCenter[0] = dotprod3(a, h);
    m_opticCenter[1] = dotprod3(a, v);
    //m_opticCenter[0] = dotprod3(a, h) - (1 + cols/2);
    //m_opticCenter[1] = (1 + rows/2) - dotprod3(a, v);
    m_fovY = 2.0 * atan2(rows, 2.0*vertical_scale);
    m_fovX = 2.0 * atan2(cols, 2.0*horizontal_scale);
    m_fx=horizontal_scale;
    m_fy=vertical_scale;
    
    // Image plane x-axis: x_axis = (H - optical_center[0]*A)/horizontal_scale;
    scalar_mult(-m_opticCenter[0], a, x_axis);
    add_vec3(h, x_axis, x_axis);
    scalar_mult(1.0/horizontal_scale, x_axis, x_axis);
    // Image plane y-axis: y_axis = (V - optical_center[1]*A)/vertical_scale;
    scalar_mult(-m_opticCenter[1], a, y_axis);
    add_vec3(v, y_axis, y_axis);
    scalar_mult(1.0/vertical_scale, y_axis, y_axis);
    
    // Set the rotation matrix R = {{H'},{-V'},{-A}}
    // Currently using R = {{H'},{V'},{A}}
    // TODO: figure out why the paper expects V' and A to be negative
    m_rotation << x_axis[0], x_axis[1], x_axis[2],
                  y_axis[0], y_axis[1], y_axis[2],
                  a[0],      a[1],      a[2];
    
    // Set the cols and rows
    m_cols = cols;
    m_rows = rows;
    
    // Copy the position
    m_position.resize(3);
    
    // Set the position of the camera
    for (int  i = 0; i < 3; i++){
	m_position[i] = c[i];
    }
    
    m_camIntrinsics << m_fx, 0, m_opticCenter[0],
      0, m_fy, m_opticCenter[1],
      0, 0, 1;
  }

  Pinhole::Pinhole(double c[3], double a[3], double h[3], double v[3], int cols, int rows, double pixelSize)
  {
        double x_axis[3] = {0.0, 0.0, 0.0}, y_axis[3] = {0.0, 0.0, 0.0};
        double AxH[3], AxV[3];
        crossprod3(a, h, AxH);
        crossprod3(a, v, AxV);
        double horizontal_scale = magnitude(AxH);
        double vertical_scale = magnitude(AxV);
        m_opticCenter[0] = dotprod3(a, h) * pixelSize;
        m_opticCenter[1] = dotprod3(a, v) * pixelSize;
        m_fovY = 2.0 * atan2(rows, 2.0*vertical_scale);
        m_fovX = 2.0 * atan2(cols, 2.0*horizontal_scale);
	m_fx=horizontal_scale * pixelSize;
        m_fy=vertical_scale * pixelSize;

        // Image plane x-axis: x_axis = (H - optical_center[0]*A)/horizontal_scale;
        scalar_mult(-m_opticCenter[0], a, x_axis);
        add_vec3(h, x_axis, x_axis);
        scalar_mult(1.0/horizontal_scale, x_axis, x_axis);
        // Image plane y-axis: y_axis = (H - optical_center[1]*A)/vertical_scale;
        scalar_mult(-m_opticCenter[1], a, y_axis);
        add_vec3(v, y_axis, y_axis);
        scalar_mult(1.0/vertical_scale, y_axis, y_axis);

        // Set the rotation matrix R = {{H'},{-V'},{-A}}
	// Currently using R = {{H'},{V'},{A}}
	// TODO: figure out why the paper expects V' and A to be negative
	m_rotation << x_axis[0], x_axis[1], x_axis[2],
		      y_axis[0], y_axis[1], y_axis[2],
		      a[0], a[1], a[2];

        // Set the cols and rows
        m_cols = cols;
        m_rows = rows;

        // Copy the position
        m_position.resize(3);

        // Set the position of the camera
        for(int  i = 0; i < 3; i++){
          m_position[i] = c[i];
        }

	m_camIntrinsics << m_fx, 0,    m_opticCenter[0],
                           0,    m_fy, m_opticCenter[1],
                           0,    0,    1;
    }

    void Pinhole::setNew(int cols, int rows, double opticCenter[2], double fx, double fy,
			 double rotation[9], double position[3])
    {
      m_position.resize(3);

      // Set the cols and rows
      m_cols = cols;
      m_rows = rows;

      // Set the optic center
      m_opticCenter[0]=opticCenter[0];
      m_opticCenter[1]=opticCenter[1];

      // Set the focal lengths
      m_fx = fx;
      m_fy = fy;

      // Set the rotation
      m_rotation<<rotation[0], rotation[1], rotation[2],
                  rotation[3], rotation[4], rotation[5],
                  rotation[6], rotation[7], rotation[8];

      // Set the position
      m_position[0] = position[0];
      m_position[1] = position[1];
      m_position[2] = position[2];

      // Set the intrinsics
      m_camIntrinsics << m_fx, 0, m_opticCenter[0],
                       0, m_fy, m_opticCenter[1],
                       0, 0, 1;
  };

  std::vector<float> Pinhole::vectorToPixel(std::vector<float> v)
  {
    std::vector<float> pixel(2);
    std::vector<float> temp(3);

    // Convert vectors to eigen matrices
    Eigen::Matrix<float, 3, 1> vec_mat;
    vec_mat << v[0], v[1], v[2];

    Eigen::Matrix<float, 3, 1> position_mat;
    position_mat << m_position[0], m_position[1], m_position[2];

    Eigen::Matrix<float, 3, 1> pix_mat;

    // Equation
    pix_mat = m_camIntrinsics * m_rotation * vec_mat - m_camIntrinsics * m_rotation * position_mat;

    // Convert eigen matrix to vector
    Eigen::Map< Eigen::Matrix<float, 3, 1> >(&temp[0], pix_mat.rows(), pix_mat.cols()) = pix_mat;
    pixel[0] = temp[0];
    pixel[1] = temp[1];

    return pixel;
  }

  std::vector<float> Pinhole::pixelToVector(std::vector<int> pixel, double lambda)
  {
	// Vector to store the result in
	std::vector<float> v;
	v.resize(3);

	// Temporary rotation matrix
	double R_temp[3][3];

	// Use equation 2.19 from http://www.epixea.com/research/multi-view-coding-thesisse8.html#x13-33005r9

	Eigen::Matrix<float, 3, 1> result;
	Eigen::Matrix<float, 3, 1> temp_result;
	Eigen::Matrix<float, 3, 1> homogeneousPixel;

	homogeneousPixel << pixel[0], pixel[1], 1;
	Eigen::Matrix<float, 3, 1> temp_position;

	temp_position[0] = m_position[0];
	temp_position[1] = m_position[1];
	temp_position[2] = m_position[2];

	v[0] = result[0];
	v[1] = result[1];
	v[2] = result[2];

	temp_result = m_rotation.inverse() * m_camIntrinsics.inverse() * homogeneousPixel;

	result = temp_position + lambda * temp_result.normalized();
	/*
	cout<<"result="<<result[0]<<", "<<result[1]<<", "<<result[2]<<endl;
	
	Eigen::Map< Eigen::Matrix<float, 3, 1> >(&v[0], result.rows(), result.cols()) = result;

	std::vector <float> finalResult;
	finalResult.resize(3);

	finalResult[0] = v[0];
	finalResult[1] = v[1];
	finalResult[2] = v[2];
	cout<<"finalResult="<<finalResult[0]<<", "<<finalResult[1]<<", "<<finalResult[2]<<endl;
	//return finalResult;
	*/

	std::vector <float> finalVec;
	finalVec.resize(3);
	finalVec[0]=result[0];
	finalVec[1]=result[1];
	finalVec[2]=result[2];

	return finalVec;
  }

  std::vector<float> Pinhole::pixelToVector(std::vector<int> pixel)
  {
	// Vector to store the result in
	std::vector<float> v;
	v.resize(3);

	// Temporary rotation matrix
	double R_temp[3][3];

	// Use equation 2.19 from http://www.epixea.com/research/multi-view-coding-thesisse8.html#x13-33005r9
	//float temp[3];
	Eigen::Matrix<float, 3, 1> result;
	Eigen::Matrix<float, 3, 1> temp_result;
	Eigen::Matrix<float, 3, 1> homogeneousPixel;

	homogeneousPixel << pixel[0], pixel[1], 1;
	Eigen::Matrix<float, 3, 1> temp_position;

	temp_position[0] = m_position[0];
	temp_position[1] = m_position[1];
	temp_position[2] = m_position[2];

	v[0] = result[0];
	v[1] = result[1];
	v[2] = result[2];

	temp_result = m_rotation.inverse() * m_camIntrinsics.inverse() * homogeneousPixel;

	result = temp_result.normalized();

	Eigen::Map< Eigen::Matrix<float, 3, 1> >(&v[0], result.rows(), result.cols()) = result;

	std::vector <float> finalResult;
	finalResult.resize(3);

	finalResult[0] = v[0];
	finalResult[1] = v[1];
	finalResult[2] = v[2];
	return finalResult;
  }

  // Writes a file with the camera parameters that can be loaded later
  // Will overwrite any file of the same name()
  void Pinhole::writeFile(std::string outputFilename)
  {
    // Open the file
    ofstream output(outputFilename.c_str());
    
    // Write all of the parameters
    output << "CAMERA_MATRIX " << m_camIntrinsics(0,0) << " " << m_camIntrinsics(0,1) << " " << m_camIntrinsics(0,2) << " " 
                               << m_camIntrinsics(1,0) << " " << m_camIntrinsics(1,1) << " " << m_camIntrinsics(1,2) << " " 
                               << m_camIntrinsics(2,0) << " " << m_camIntrinsics(2,1) << " " << m_camIntrinsics(2,2) << endl;

    output << "ROTATION " << m_rotation(0,0) << " " << m_rotation(0,1) << " " << m_rotation(0,2) << " " 
                          << m_rotation(1,0) << " " << m_rotation(1,1) << " " << m_rotation(1,2) << " " 
                          << m_rotation(2,0) << " " << m_rotation(2,1) << " " << m_rotation(2,2) << endl;

    output << std::setprecision(12) << "TRANSLATION " << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;

    output << "WIDTH_HEIGHT " << m_cols << " " << m_rows << endl;

    if(m_distortion.size() != 0)
    {
      output << "DISTORTION_COEFFICIENTS";
      for(std::vector<float>::iterator it = m_distortion.begin(); it != m_distortion.end(); it++){
        output << " " << *it;
      }
      output << endl;
    }
  }

  void Pinhole::print()
  {

    // Write all of the parameters
    cout<< "CAMERA_MATRIX " << m_camIntrinsics(0,0) << " " << m_camIntrinsics(0,1) << " " << m_camIntrinsics(0,2) << " " 
	<< m_camIntrinsics(1,0) << " " << m_camIntrinsics(1,1) << " " << m_camIntrinsics(1,2) << " " 
	<< m_camIntrinsics(2,0) << " " << m_camIntrinsics(2,1) << " " << m_camIntrinsics(2,2) << endl;
    
    cout << "ROTATION " << endl;
    cout << m_rotation(0,0) << " " << m_rotation(0,1) << " " << m_rotation(0,2) << endl;
    cout << m_rotation(1,0) << " " << m_rotation(1,1) << " " << m_rotation(1,2) << endl; 
    cout << m_rotation(2,0) << " " << m_rotation(2,1) << " " << m_rotation(2,2) << endl;
    
    cout<< "TRANSLATION " << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;

    cout << "WIDTH_HEIGHT " << m_cols << " " << m_rows << endl;

    if(m_distortion.size() != 0){
      cout << "DISTORTION_COEFFICIENTS";
      for(std::vector<float>::iterator it = m_distortion.begin(); it != m_distortion.end(); it++){
        cout << " " << *it;
      }
      cout << endl;
    }
  }

  void Pinhole::GetRotation (double get_rotation[3][3])
  {
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++) 
	//Eigen::Matrix structures are column major, so they need to be transposed
	get_rotation[j][i] = m_rotation(3*i + j);
  };
  
  void  Pinhole::GetIntrinsics(double get_intrinsics[3][3])
  {
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++)
	//Eigen::Matrix structures are column major, so they need to be transposed
	get_intrinsics[j][i] = m_camIntrinsics(3*i + j);
  };
  
  void  Pinhole::GetModel(double get_optical_center[2], double* get_fov_x, double* get_fov_y)
  {
    // Set the optical center
    get_optical_center[0] = m_opticCenter[0];
    get_optical_center[1] = m_opticCenter[1];
    
    // Set the fov
    *get_fov_x = m_fovX;
    *get_fov_y = m_fovY;
  }
    
  void  Pinhole::GetSize(int& get_cols, int& get_rows)
  {
    get_cols = m_cols;
    get_rows = m_rows;
  }
  
  void  Pinhole::GetPosition(std::vector<float>& get_position)
  {
    get_position.resize(3);
    
    for (int i = 0; i < 3; i++){
      get_position[i] = m_position[i];
    }
  }
  
  bool  Pinhole::GetDistortion(std::vector<float>& get_distortion)
  {
    get_distortion.clear();
    
    if(m_distortion.size() == 0){
      return false;
    }
    
    else{
      for(std::vector<float>::iterator it = m_distortion.begin(); it != m_distortion.end(); it++){
	get_distortion.push_back(*it);
      }
      return true;
    }
  }
  
  // mosaic_processing needs this function.
  void Pinhole::
  GetPinholeModel(double C[3], double A[3], double H[3], double V[3],
                  int cols, int rows, float pixelSize, double optical_center[2],
		  double *fov_X, double *fov_Y)
    {
        double x_axis[3] = {0.0, 0.0, 0.0}, y_axis[3] = {0.0, 0.0, 0.0};
        double AxH[3], AxV[3];//, optical_center[2];
        crossprod3(A, H, AxH);
        crossprod3(A, V, AxV);
        double horizontal_scale = magnitude(AxH);
        double vertical_scale = magnitude(AxV);
        optical_center[0] = dotprod3(A, H);
        optical_center[1] = dotprod3(A, V);
        /*double*/ *fov_Y = 2.0 * atan2(rows, 2.0*vertical_scale);
        /*double*/ *fov_X = 2.0 * atan2(cols, 2.0*horizontal_scale);
        /*
          cout << "Horizontal scale = " << horizontal_scale << endl;
          cout << "Vertical scale = " << vertical_scale << endl;

          cout << "Focal length = " << ((horizontal_scale + vertical_scale)*0.5 *
          pixelSize)
          << endl;

          cout << "Optical center pixel = [" << lround(optical_center[0]) << ", "
          << lround(optical_center[1]) << "]" << endl;
        */
        // Image plane x-axis: x_axis = (H - optical_center[0]*A)/horizontal_scale;
        scalar_mult(-optical_center[0], A, x_axis);
        add_vec3(H, x_axis, x_axis);
        scalar_mult(1.0/horizontal_scale, x_axis, x_axis);
        // Image plane y-axis: y_axis = (H - optical_center[1]*A)/vertical_scale;
        scalar_mult(-optical_center[1], A, y_axis);
        add_vec3(V, y_axis, y_axis);
        scalar_mult(1.0/vertical_scale, y_axis, y_axis);
        /*
          cout << "Center of projection = {"
          << C[0] << ", "
          << C[1] << ", "
          << C[2] << "}" << endl;

          cout << "Image plane X-axis = {"
          << x_axis[0] << ", "
          << x_axis[1] << ", "
          << x_axis[2] << "}" << endl;

          cout << "Image plane Y-axis = {"
          << y_axis[0] << ", "
          << y_axis[1] << ", "
          << y_axis[2] << "}" << endl;

          cout << "Image plane Z(A)-axis = {"
          << A[0] << ", "
          << A[1] << ", "
          << A[2] << "}" << endl;

          cout << "FOV X (deg) = " << *fov_X * 180.0 / M_PI << endl;
          cout << "FOV Y (deg) = " << *fov_Y * 180.0 / M_PI << endl;
        */
    }
}
