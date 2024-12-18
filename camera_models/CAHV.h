#ifndef _CAHV_H_
#define _CAHV_H_

#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
using namespace std;

namespace at 
{

    class CAHV
    {
    public:

        long int m_height;
        long int m_width;
        Eigen::Vector3f m_c;
        Eigen::Vector3f m_a;
        Eigen::Vector3f m_h;
        Eigen::Vector3f m_v;
        vector<float> m_offset;
        vector<float> m_quaternion;

        CAHV();
        CAHV(vector<float> c, vector<float> a, vector<float> h, vector<float> v,
             vector<float>offset, long int width, long int height);
        CAHV(vector<float> c, vector<float> a, vector<float> h, vector<float> v,
             vector<float> offset, vector<float> quaternion,
             long int width, long int height);
	CAHV(Eigen::Vector2d focalLength, Eigen::Vector2d opticalCenter, int imgWidth, int imgHeight, 
	     Eigen::Matrix3d rotation, Eigen::Vector3d translation);

	// This constructor will throw exceptions when loading fails
        // This will delete the object
        // Will apply a quaternion if one is found
        // 0 means that the file was not found
        // 1 means that the file was misconfigured
	CAHV(string cahv_calibration);
        ~CAHV() { };
    
        vector<float> pixelToVector(vector<int> thisPixel);
        vector<int>   vectorToPixel(vector<float> thisVector);

	void getImagePlane(double a_unit[3], double h_unit[3], double v_unit[3]);

	// Writes the current CAHV model to a file so that it can be loaded later
        // Saves final CAHV, will not save a quaternion if one was used to create this object
	void writeFile(std::string outputFilename);
	void print();

    private:

	// From: http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
	// These functions split a string given a delimiter
	// Should these go to common/string_util
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	    std::stringstream ss(s);
	    std::string item;
	    while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	    }
	    return elems;
	}


	std::vector<std::string> split(const std::string &s, char delim = ' ') {
	    std::vector<std::string> elems;
	    split(s, delim, elems);
	    return elems;
	}
    };
}

#endif
