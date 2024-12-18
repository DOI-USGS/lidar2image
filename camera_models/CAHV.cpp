#include <iostream>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <iomanip>

#include "CAHV.h"
#include "../common/GeometricTransf.h"

using namespace std;

namespace at
{
    CAHV::CAHV()
    {

    }

    //not sure if this constructor is needed.
    CAHV::CAHV(vector<float> c, vector<float> a, vector<float> h, vector<float> v,
	     vector<float> offset, long int width, long int height)
    {
      m_height = height;
      m_width = width;
      m_offset.resize(3);
      m_c.resize(3);
      m_a.resize(3);
      m_h.resize(3);
      m_v.resize(3);
      for (int i=0; i < 3; i++){
	m_c[i]=c[i]+offset[i];
	m_a[i]=a[i];
	m_h[i]=h[i];
	m_v[i]=v[i];
	m_offset[i]=offset[i];
      }
    }

    CAHV::CAHV(vector<float> c, vector<float> a, vector<float> h, vector<float> v,
               vector<float> offset, vector<float> quaternion,
               long int width, long int height)
    {
        m_height = height;
        m_width = width;
        m_offset.resize(3);
        m_quaternion.resize(4);
        m_c.resize(3);
        m_a.resize(3);
        m_h.resize(3);
        m_v.resize(3);

        double c_in[3];
        double a_in[3];
        double h_in[3];
        double v_in[3];
        double temp_quaternion[4];

        for (int i=0; i < 3; i++){
	  c_in[i]=c[i];
	  a_in[i]=a[i];
	  h_in[i]=h[i];
	  v_in[i]=v[i];
	  m_offset[i]=offset[i];
        }
	
        for (int i=0; i < 4; i++){
            m_quaternion[i] = quaternion[i];
            temp_quaternion[i] = quaternion[i];
        }
  
        double c_out[3]; 
        rotateWithQuaternion( c_in, c_out, temp_quaternion);
        double a_out[3]; 
        rotateWithQuaternion( a_in, a_out, temp_quaternion);
        double h_out[3];
        rotateWithQuaternion( h_in, h_out, temp_quaternion);
        double v_out[3];
        rotateWithQuaternion( v_in, v_out, temp_quaternion);

        for (int i=0; i < 3; i++){
	  m_c[i]=c_out[i] + m_offset[i];//Sept 02, 2014 with Larry
	  m_a[i]=a_out[i];
	  m_h[i]=h_out[i];
	  m_v[i]=v_out[i];	
        }
    }

    CAHV::CAHV(string cahv_calibration)
    {
	std::ifstream fin (cahv_calibration.c_str());
	std::string line;

	bool found_offset = false;
	bool found_quaternion = false;
	int found = 0;

	double c[3],  a[3], h[3], v[3];
	m_quaternion.resize(4);
	m_offset.resize(3);

	m_height = -1;
	m_width = -1;

	// Check for proper loading
	if (fin.is_open()){// File loaded correctly
	  // Continue loading until the end of the file
	  while (!fin.eof()){
	    line.clear();
	    std::getline(fin, line);
	    
	    if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines

	      std::stringstream lineStream(line);
	      std::string keyword("");
	      lineStream >> keyword;
	      
	      if (keyword.compare(std::string("WIDTH_HEIGHT"))==0){
		found++;
	       
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 3) {// Not enough parameters
		  // Throw failed load
		  throw 1;
		}
		
		// Get the parameters
		m_width = atof(contents[1].c_str());
		m_height = atof(contents[2].c_str());
	      }
	      else if (keyword.compare(std::string("C"))==0){
		  found++;
		  
		  // Check number of parameters
		  vector<string> contents = split(line);
		  if (contents.size() < 4){// Not enough parameters
		    // Throw failed load
		    throw 1;
		  }
		  
		  // Get the parameters
		  c[0] = atof(contents[1].c_str());
		  c[1] = atof(contents[2].c_str());
		  c[2] = atof(contents[3].c_str());	
	      }
	      else if (keyword.compare(std::string("A"))==0){
		found++;
		lineStream>> a[0] >> a[1] >> a[2];
		
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 4){ // Not enough parameters
		    // Throw failed load
		    throw 1;
		  }
		
		// Get the parameters
		a[0] = atof(contents[1].c_str());
		a[1] = atof(contents[2].c_str());
		a[2] = atof(contents[3].c_str());
	      }
	      else if (keyword.compare(std::string("H"))==0){
		found++;
		lineStream >> h[0] >> h[1] >> h[2];
		
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 4){// Not enough parameters
		  // Throw failed load
		  throw 1;
		}
		
		// Get the parameters
		h[0] = atof(contents[1].c_str());
		h[1] = atof(contents[2].c_str());
		h[2] = atof(contents[3].c_str());
	      }
	      else if (keyword.compare(std::string("V"))==0){
		found++;		  
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 4){// Not enough parameters 
		  // Throw failed load
		  throw 1;
		}
		
		// Get the parameters
		v[0] = atof(contents[1].c_str());
		v[1] = atof(contents[2].c_str());
		v[2] = atof(contents[3].c_str());
	      }
	      else if (keyword.compare(std::string("QUATERNION"))==0){
		found_quaternion = true;	  
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 5){ // Not enough parameters
		  // Throw failed load
		  throw 1;
		}
		
		// Get the parameters
		m_quaternion[0] = atof(contents[1].c_str());
		m_quaternion[1] = atof(contents[2].c_str());
		m_quaternion[2] = atof(contents[3].c_str());
		m_quaternion[3] = atof(contents[4].c_str());
	      }
	      else if (keyword.compare(std::string("OFFSET"))==0){

		found_offset = true;
		
		// Check number of parameters
		vector<string> contents = split(line);
		if (contents.size() < 4){// Not enough parameters
		  // Throw failed load
		  throw 1;
		}
		
		// Get the parameters
		m_offset[0] = atof(contents[1].c_str());
		m_offset[1] = atof(contents[2].c_str());
		m_offset[2] = atof(contents[3].c_str());
	      }
	    }
	  }
	  
	  // Check for correct number of parameters
	  if (found < 5){// Not enough parameters
	    // Throw failed load exception
	    throw 1;
	  }
	  
	  // Add offset to the current position if an offset was loaded
	  if (found_offset){// Offset was found
	    // Add the offset to the current position
	    for (int i = 0; i < 3; i++){
	      m_c[i] = c[i] + m_offset[i];
	    }
	  }
	  
	  else{ // No offset found
	    // Use c without modification
	    for (int i = 0; i < 3; i++){
	      m_c[i] = c[i];
	    }
	  }
	  
	  // Apply quaternion if one was loaded
	  if (found_quaternion){ // Quaternion was found
	    // Copy the quaternion
	    double temp_quaternion[4];
	    for (int i = 0; i < 4; i++){
	      temp_quaternion[i] = m_quaternion[i];
	    } 
	    
	    // Apply the quaternion
	    double c_out[3]; 
	    rotateWithQuaternion( c, c_out, temp_quaternion);
	    double a_out[3]; 
	    rotateWithQuaternion( a, a_out, temp_quaternion);
	    double h_out[3];
	    rotateWithQuaternion( h, h_out, temp_quaternion);
	    double v_out[3];
	    rotateWithQuaternion( v, v_out, temp_quaternion);
	    
	    // Save the results as the new cavh parameters
	    for(int i = 0; i < 3; i++){
	      m_c[i] = c_out[i];
	      m_a[i] = a_out[i];
	      m_h[i] = h_out[i];
	      m_v[i] = v_out[i];
	    }
	  }
	  
	  else{ // No quaternion was found
	    // Use parameters without modification
	    for(int i = 0; i < 3; i++){
	      m_c[i] = c[i];
	      m_a[i] = a[i];
	      m_h[i] = h[i];
	      m_v[i] = v[i];
	    }
	  }
	}
	
	else{// File did not load correctly
	  // Throw failed open exception
	  throw 0;
	}
    }
  
    //build cahv (in pixels) from pinhole
    //focalLength, opticalCenter, imgWidth and imgHeight in pixels
    CAHV::CAHV(Eigen::Vector2d focalLength, Eigen::Vector2d opticalCenter, int imgWidth, int imgHeight, 
	       Eigen::Matrix3d rotation, Eigen::Vector3d translation)
    {
      
      //m_c = rotation*centerProj + translation;
      for (int i=0; i < 3; i++){  
	m_c[i]=translation[i];  
      }

      for (int i = 0; i < 3; i++){
	m_a[i] = rotation(2,i);
      }

      //vr_real32 fH = f/pixelSize[0], fV = f/pixelSize[1];

      for (int i = 0; i < 3; i++){
	m_h[i] = focalLength[0] * rotation(0,i) + opticalCenter[0] * m_a[i];
      }
      for (int i = 0; i < 3; i++){
	m_v[i] = focalLength[1] * rotation(1,i) + opticalCenter[1] * m_a[i];
      }

      m_height = (long int)imgHeight;
      m_width = (long int)imgWidth;
    }

    vector<float> CAHV::pixelToVector(vector<int> thisPixel)
    {
        //vertical component
        Eigen::Vector3f va = m_v + (-thisPixel[1])*m_a;
        //horizontal component
        Eigen::Vector3f vb = m_h + (-thisPixel[0])*m_a;
        //cross product
        Eigen:: Vector3f pVec = va.cross(vb);
        //normalization
        //cout<<"pVec="<<pVec<<endl;
        //cout<<"sqnorm="<<sqrt(pVec.squaredNorm())<<endl;
        pVec = pVec / pVec.norm();
        //cout<<"nVec="<<pVec<<endl;

        // The vector should be pointing in the same directions as A, if it
        // isn't, flip it:
        if (pVec.dot(m_a) < 0.0)
            pVec *= -1;

        vector<float> thisVector(3);

        thisVector[0] = pVec[0];
        thisVector[1] = pVec[1];
        thisVector[2] = pVec[2];


        return thisVector;
    }


  vector<int>
  CAHV::vectorToPixel(vector<float> thisVector)
  {
    vector<int> thisPixel(2);
    Eigen::Vector3f pVec(thisVector.data());
    double dDot = pVec.dot(m_a);

    // Ara probably wants to make this condition > 0.0... just sayin'
    // -- LJE
    if (dDot != 0.0)
    {
      thisPixel[0] = pVec.dot(m_h) / dDot;
      thisPixel[1] = pVec.dot(m_v) / dDot;
    }
    else
    {
      thisPixel[0] = numeric_limits<int>::max();
      thisPixel[1] = numeric_limits<int>::max();
    }

    return thisPixel;
  }

  void CAHV::getImagePlane(double a_unit[3], double h_unit[3], double v_unit[3])
  {
    // Extract image plane vectors from H and V
    Eigen::Vector3f h_img_plane = (m_h - m_a.dot(m_h) * m_a) / (m_a.cross(m_h).norm());
    Eigen::Vector3f v_img_plane = (m_v - m_a.dot(m_v) * m_a) / (m_a.cross(m_v).norm());

    // Copy to the output vectors
    for(int i = 0; i < 3; i++)
    {
      a_unit[i] = m_a(i);
      h_unit[i] = h_img_plane(i);
      v_unit[i] = v_img_plane(i);
    }
  }

  // Writes the current CAHV model to a file so that it can be loaded later
  void CAHV::writeFile(string outputFilename)
  {
    // Open the file
    ofstream output(outputFilename.c_str());

    // Write all of the parameters
    output << std::setprecision(8) << "C " << m_c[0] << " " << m_c[1] << " " << m_c[2] << endl;
    output << std::setprecision(8) << "A " << m_a[0] << " " << m_a[1] << " " << m_a[2] << endl;
    output << std::setprecision(8) << "H " << m_h[0] << " " << m_h[1] << " " << m_h[2] << endl;
    output << std::setprecision(8) << "V " << m_v[0] << " " << m_v[1] << " " << m_v[2] << endl;

    output << "WIDTH_HEIGHT " << m_width << " " << m_height << endl;
  }
  void CAHV::print()
  {
    cout << "CAHV model params: " << endl;
    cout << "c = " << m_c[0] << ", " << m_c[1] << ", " << m_c[2] << endl;
    cout << "a = " << m_a[0] << ", " << m_a[1] << ", " << m_a[2] << endl;
    cout << "h = " << m_h[0] << ", " << m_h[1] << ", " << m_h[2] << endl;
    cout << "v = " << m_v[0] << ", " << m_v[1] << ", " << m_v[2] << endl;
  }

}
