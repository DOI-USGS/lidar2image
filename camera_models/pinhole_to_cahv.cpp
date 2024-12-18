#include <vector>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <iomanip>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


#include "CAHV.h"
#include "pinhole.h"

using namespace std;
using namespace at;

int main( int argc, char *argv[] )
{
  if (argc<3){
    cout<<"Usage: pinhole_to_cahv cameraFilename resultsDirname"<<endl;
    exit(1);
  }

   string cameraFilename = string(argv[1]);
   string resultsDirname = string(argv[2]);
   
   cout << "Attempting to load file as pinhole model first" << endl;
   Pinhole pin(cameraFilename);
   cout << "Loaded pinhole model" << endl;
 
   //convert pinhole to CAHV
   cout<<"Converting to CAHV ..."<<endl;

   double initTranslation[3];
   initTranslation[0] = pin.GetCameraPositionX();
   initTranslation[1] = pin.GetCameraPositionY();
   initTranslation[2] = pin.GetCameraPositionZ();
   cout<<"initTranslation="<<initTranslation[0]<<", "
       <<initTranslation[1]<<", "<<initTranslation[2]<<endl;
   
   double input_rotation[3][3];
   pin.GetRotation(input_rotation);
   cout<<"input rotation="<<endl;
   cout<<input_rotation[0][0]<<", "<<input_rotation[0][1]<<", "<<input_rotation[0][2]<<endl;
   cout<<input_rotation[1][0]<<", "<<input_rotation[1][1]<<", "<<input_rotation[1][2]<<endl;
   cout<<input_rotation[2][0]<<", "<<input_rotation[2][1]<<", "<<input_rotation[2][2]<<endl;
   
   Eigen::Vector2d focalLength;
   focalLength[0] = (float)(pin.GetFocalLengthHor());
   focalLength[1] = (float)(pin.GetFocalLengthVer());
   cout<<"focalLength (in pixels)="<<focalLength[0]<<", "<<focalLength[1]<<endl;

   Eigen::Vector2d opticalCenter;
   opticalCenter[0] = (float)(pin.GetOpticCenterHor());
   opticalCenter[1] = (float)(pin.GetOpticCenterVer());
   cout<<"opticalCenter (in pixels)="<<opticalCenter[0]<<", "<<opticalCenter[1]<<endl;

   Eigen::Matrix3d rotation;
   rotation(0,0)=input_rotation[0][0]; rotation(0,1)=input_rotation[0][1]; rotation(0,2)=input_rotation[0][2];
   rotation(1,0)=input_rotation[1][0]; rotation(1,1)=input_rotation[1][1]; rotation(1,2)=input_rotation[1][2];
   rotation(2,0)=input_rotation[2][0]; rotation(2,1)=input_rotation[2][1]; rotation(2,2)=input_rotation[2][2];

   Eigen::Vector3d translation;
   translation(0)=initTranslation[0];
   translation(1)=initTranslation[1];
   translation(2)=initTranslation[2];

   vector<int> pixel; pixel.resize(2);
   pixel[0]=opticalCenter[0]+200;
   pixel[1]=opticalCenter[1]+200;

   vector<float> vector3d_pinhole = pin.pixelToVector(pixel, 1);
   cout<<"Pinhole Vector3D="<<vector3d_pinhole[0]<<", "<<vector3d_pinhole[1]<<", "<<vector3d_pinhole[2]<<endl;
   
   CAHV outCAHV = CAHV(focalLength, opticalCenter, (long int)pin.GetCols(), (long int)pin.GetRows(), 
		       rotation, translation);
   
   cout << "CAHV camera model out: " << endl;
   cout << "\tf_c = " << outCAHV.m_c[0] << ", " << outCAHV.m_c[1] << ", " << outCAHV.m_c[2] << endl;
   cout << "\tf_a = " << outCAHV.m_a[0] << ", " << outCAHV.m_a[1] << ", " << outCAHV.m_a[2] << endl;
   cout << "\tf_h = " << outCAHV.m_h[0] << ", " << outCAHV.m_h[1] << ", " << outCAHV.m_h[2] << endl;
   cout << "\tf_v = " << outCAHV.m_v[0] << ", " << outCAHV.m_v[1] << ", " << outCAHV.m_v[2] << endl;
   cout<< "\tf_width height = "<<outCAHV.m_width<<", "<<outCAHV.m_height<<endl;

   vector<float> vector3d_cahv=outCAHV.pixelToVector(pixel);
   cout<<"CAHV Vector3D="<<outCAHV.m_c[0]+vector3d_cahv[0]
       <<", "<<outCAHV.m_c[1]+vector3d_cahv[1]
       <<", "<<outCAHV.m_c[2]+vector3d_cahv[2]<<endl;
 
   string resultCamFilename = resultsDirname+"/cahv.txt";   
   outCAHV.writeFile(resultCamFilename);
  
}







