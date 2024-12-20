/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d_omp.h>

#include "pc.h"
#include "Tiling.h"

using namespace std;

namespace at
{

    pcl::PointXYZRGB readLineFromPCFile(string line, int r, int g, int b)
    {
        pcl::PointXYZRGB thisPoint;
        float x = 0, y = 0, z = 0; 
        int l_r=r; 
        int l_g=g; 
        int l_b=b;
        sscanf(line.c_str(), "%f %f %f %d %d %d", &x, &y, &z, &l_r, &l_g, &l_b);
        thisPoint.x=x;
        thisPoint.y=y;
        thisPoint.z=z;
        uint32_t rgbValue = (static_cast<uint8_t>(l_r) << 16 | static_cast<uint8_t>(l_g) << 8 | static_cast<uint8_t>(l_b));
        thisPoint.rgb = *reinterpret_cast<float*>(&rgbValue);
        return thisPoint;
    }
  
    void readHeaderFromPCFile(string line, int *width, int *height, int *hasHeader)
    {
        string widthString, heightString, pcWidth, pcHeight;
        stringstream ssline(line);
        ssline>>widthString>>pcWidth>>heightString>>pcHeight;
        if ((widthString.compare("WIDTH:")==0)&&(heightString.compare("HEIGHT:")==0)){
            *hasHeader = 1;
            *height=atoi(pcHeight.c_str());
            *width=atoi(pcWidth.c_str());
        }
        else{
            *hasHeader=0;
            *height=1;
            *width=0;
        }
    }
    //read the poincloud filename with colors if they exist otherwise with fixed rgb colors.
    void appendCloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, string pointCloudFilename, int subSampleStep, int r, int g, int b)
    {
        string line;
        ifstream file (pointCloudFilename.c_str());
        pcl::PointXYZRGB thisPoint;
        int index = 0;

        int width, height, hasHeader;

        if (file.is_open()){
            //check the header
            getline (file,line);
            readHeaderFromPCFile(line, &width, &height, &hasHeader);
            if (!hasHeader){
                thisPoint=readLineFromPCFile(line, r, g, b);
                if (((thisPoint.x!=0.0)||(thisPoint.y!=0.0)||(thisPoint.z!=0.0)) && (subSampleStep*(index/subSampleStep) == index) ){
                    cloud->points.push_back(thisPoint);
                }
                index++;
            }
            //loop over the remaining lines
            while ( file.good() ){
                getline (file,line);
                thisPoint=readLineFromPCFile(line, r, g, b);
                if (((thisPoint.x!=0.0)||(thisPoint.y!=0.0)||(thisPoint.z!=0.0)) && (subSampleStep*(index/subSampleStep) == index) ){
                    cloud->points.push_back(thisPoint);
                }
                index++;
            }
            cout<<"numPoints="<<cloud->points.size()<<endl;
            file.close();
            if (hasHeader==1){
                cloud->width=width;
                cloud->height=height;
            }
            if (hasHeader==0){
                cloud->height=1;
                cloud->width=cloud->points.size();
            }
        }
        else{
            cout << "Unable to open file." << endl;
        } 
        return;
    }

    //read the poincloud filename
    void appendCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, string pointCloudFilename, int subSampleStep)
    {
        string line;
        ifstream file (pointCloudFilename.c_str());
        pcl::PointXYZ thisPoint;
        int index = 0;
        cout<<"numPointsBefore="<<cloud->points.size()<<endl;
        if (file.is_open()){
            while ( file.good() ){
                getline (file,line);
                float x = 0, y = 0, z = 0; 
//                int r, g, b;
                //sscanf(line.c_str(), "%f %f %f %d %d %d", &x, &y, &z, &r, &g, &b);
                sscanf(line.c_str(), "%f %f %f", &x, &y, &z);
                if (((x!=0.0)||(y!=0.0)||(z!=0.0)) && (subSampleStep*(index/subSampleStep) == index) ){	
                    thisPoint.x=x;
                    thisPoint.y=y;
                    thisPoint.z=z;
                    cloud->points.push_back(thisPoint);
                    cloud->width = cloud->width+1;
                }
                index++;
            }
            cout<<"append_numPoints="<<cloud->points.size()<<endl;
            file.close();
        }
        else{
            cout << "Unable to open file." << endl;
        } 
        return;
    }


    // functions below this line will be removed
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr  ReadPCFile(string pointCloudFilename, int numPoints, int subSampleStep)
    {
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
        cloud->width = numPoints/subSampleStep;
        cloud->height = 1;
        //cloud->is_dense = false;
        cloud->points.resize (cloud->width * cloud->height);
        int counter = 0;
        ifstream outFile (pointCloudFilename.c_str());
        string line;
        if (outFile.is_open())
        {
            while ( outFile.good() )
            {
                getline (outFile,line);
                float x = 0, y = 0, z = 0, r = 255, g = 0, b = 0;
                sscanf(line.c_str(), "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
                if ((x!=0.0)||(y!=0.0)||(z!=0.0))
                {
                    cloud->points[counter].x = x;
                    cloud->points[counter].y = y;
                    cloud->points[counter].z = z;
                    cloud->points[counter].r = r;
                    cloud->points[counter].g = g;
                    cloud->points[counter].b = b;
	      
                    if ((r<10) && (g<10) && (b<10)){
                        //if points are too dark show in red to distinguish from background
                        cloud->points[counter].r = 255;
                        cloud->points[counter].g = 0;
                        cloud->points[counter].b = 0;
                    }
	      
                    counter++;
                }
                for(int i=0; i<subSampleStep-1; i++)
                    getline(outFile, line);
            }
            outFile.close(); 
        }
        return cloud;
    }

    void saveOrganizedPointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, string filename)
    {
      cout<<"filename="<<filename<<endl;
      ofstream file(filename.c_str());
      file<<"WIDTH: "<<pointCloud->width<<" HEIGHT:"<<pointCloud->height<<endl;
      for (unsigned int i = 0; i < pointCloud->points.size(); i++){
	int r = static_cast<int>(pointCloud->points[i].r);
	int g = static_cast<int>(pointCloud->points[i].g);
	int b = static_cast<int>(pointCloud->points[i].b);
	file<<std::fixed /*<< std::setprecision(10)*/
	    <<pointCloud->points[i].x<<" "<<pointCloud->points[i].y<<" "<<pointCloud->points[i].z<<" "
	    <<r<<" "<<g<<" "<<b<<endl;
      }
      file.close();
    }
  
    void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, string filename)
    {
      cout<<"filename="<<filename<<endl;
      ofstream file(filename.c_str());
      for (unsigned int i = 0; i < pointCloud->points.size(); i++){
	int r = static_cast<int>(pointCloud->points[i].r);
	int g = static_cast<int>(pointCloud->points[i].g);
	int b = static_cast<int>(pointCloud->points[i].b);
	file<<std::fixed 
	    <<pointCloud->points[i].x<<" "<<pointCloud->points[i].y<<" "<<pointCloud->points[i].z<<" "
	    <<r<<" "<<g<<" "<<b<<endl;
      }
      file.close();
    }
  
    void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, float xOffset, float yOffset, string filename)
    {

      ofstream file(filename.c_str());
      file<<"WIDTH: "<<pointCloud->width<<" HEIGHT: "<<pointCloud->height<<endl;
      for (unsigned int i = 0; i < pointCloud->points.size(); i++){
	int r = static_cast<int>(pointCloud->points[i].r);
	int g = static_cast<int>(pointCloud->points[i].g);
	int b = static_cast<int>(pointCloud->points[i].b);
	file<<std::fixed 
	    <<pointCloud->points[i].x + xOffset<<" "<<pointCloud->points[i].y+yOffset<<" "<<pointCloud->points[i].z<<" "
	    <<r<<" "<<g<<" "<<b<<endl;
      }
      file.close();
    }
  
  void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      int r = static_cast<int>(pointCloud.points[i].r);
      int g = static_cast<int>(pointCloud.points[i].g);
      int b = static_cast<int>(pointCloud.points[i].b);
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<" "
	  <<r<<" "<<g<<" "<<b<<endl;
    }
    file.close();
  }

 

  void savePointNormal(pcl::PointCloud<pcl::PointNormal> pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      int r = 0;//static_cast<int>(pointCloud.points[i].r);
      int g = 0;//static_cast<int>(pointCloud.points[i].g);
      int b = 0;//static_cast<int>(pointCloud.points[i].b);
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<" "
	  <<r<<" "<<g<<" "<<b<<endl;
    }
    file.close();
  }

  void saveOrganizedPointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    file<<"WIDTH: "<<pointCloud.width<<" HEIGHT:"<<pointCloud.height<<endl;
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      int r = static_cast<int>(pointCloud.points[i].r);
      int g = static_cast<int>(pointCloud.points[i].g);
      int b = static_cast<int>(pointCloud.points[i].b);
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<" "
	  <<r<<" "<<g<<" "<<b<<endl;
    }
    file.close();
  }
  
  void appendPointCloudRGBToFile(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str(), ios::app);
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      int r = static_cast<int>(pointCloud.points[i].r);
      int g = static_cast<int>(pointCloud.points[i].g);
      int b = static_cast<int>(pointCloud.points[i].b);
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<" "
	  <<r<<" "<<g<<" "<<b<<endl;
    }
    file.close();
  }
  
  void savePointCloudPtr(pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    for (unsigned int i = 0; i < pointCloud->points.size(); i++){
      file<<std::fixed<<pointCloud->points[i].x<<" "<<pointCloud->points[i].y<<" "<<pointCloud->points[i].z<<endl;
    }
    file.close();
  }
  void savePointCloud(pcl::PointCloud<pcl::PointXYZ> pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<endl;
    }
    file.close();
  }
  void savePointCloud(const pcl::PointCloud<pcl::PointXYZ>& pointCloud, string filename)
  {
    cout<<"filename="<<filename<<endl;
    ofstream file(filename.c_str());
    for (unsigned int i = 0; i < pointCloud.points.size(); i++){
      file<<std::fixed<<pointCloud.points[i].x<<" "<<pointCloud.points[i].y<<" "<<pointCloud.points[i].z<<endl;
    }
    file.close();
  }
  
  pcl::PointCloud<pcl::PointXYZ>::Ptr  ReadPCNoRGBFile(string pointCloudFilename, int numPoints, int subSampleStep)
  {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    cloud->width = numPoints/subSampleStep;
    cloud->height = 1;
    //cloud->is_dense = false;
    cloud->points.resize (cloud->width * cloud->height);
    int counter = 0;
    ifstream outFile (pointCloudFilename.c_str());
    string line;
    if (outFile.is_open()){
      while ( outFile.good() ){
	getline (outFile,line);
	float x = 0, y = 0, z = 0;
	sscanf(line.c_str(), "%f %f %f", &x, &y, &z);
	if ((x!=0.0)||(y!=0.0)||(z!=0.0)){
	  cloud->points[counter].x = x;
	  cloud->points[counter].y = y;
	  cloud->points[counter].z = z; 
	  counter++;
	}
	for(int i=0; i<subSampleStep-1; i++){
	  getline(outFile, line);
	}
      }
      outFile.close(); 
    }
    return cloud;
  }
  
  
  int CountEntries(string pointCloudFilename)
  {
    string line;
    ifstream myfile (pointCloudFilename.c_str());
    
    int numPoints = 0;
    if (myfile.is_open()){
      while ( myfile.good() ){
	getline (myfile,line);
	float x = 0, y = 0, z = 0, r = 0, g = 0, b = 0;
	sscanf(line.c_str(), "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
	if ((x!=0.0)||(y!=0.0)||(z!=0.0)){
	    numPoints++;
	}
      }
      myfile.close();
    }
    
    else cout << "Unable to open file." << endl;
    cout<<"numPoints="<<numPoints<<endl;
    return numPoints;
  }
  
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr makeTransformedPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, Eigen::Matrix4f transformF)
  {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr transformPC (new pcl::PointCloud<pcl::PointXYZRGB>);

    Eigen::Matrix4d transform;
    for (int  i =0; i < 4; i++){
      for (int j=0; j < 4; j++){
	transform(i,j)=(double)transformF(i,j);
      }
    }

    Eigen::Vector4d initPoint;
    Eigen::Vector4d finalPoint;
    pcl::PointXYZRGB point;

    for (unsigned int i=0; i < pointCloud->points.size(); i++){
      
      initPoint(0) = (double)pointCloud->points[i].x;
      initPoint(1) = (double)pointCloud->points[i].y;
      initPoint(2) = (double)pointCloud->points[i].z;
      initPoint(3) = 1.0;

      finalPoint = transform*initPoint;
          
      point.x = (float)finalPoint(0);
      point.y = (float)finalPoint(1);
      point.z = (float)finalPoint(2);
      point.r = pointCloud->points[i].r;
      point.g = pointCloud->points[i].g;
      point.b = pointCloud->points[i].b;

      transformPC->points.push_back(point);    
    }
    transformPC->width = pointCloud->width;
    transformPC->height = pointCloud->height;

    return transformPC;  
  }


  struct BoundingBox  GetBoundingBox(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, pcl::PointXYZRGB noDataValue)
  { 

    BoundingBox bbox3;
    bbox3.first.resize(3);
    bbox3.last.resize(3);

    bbox3.first[0] = 10000000.0;
    bbox3.first[1] = 10000000.0;
    bbox3.first[2] = 10000000.0;
    bbox3.last[0]  =-10000000.0;
    bbox3.last[1]  =-10000000.0;
    bbox3.last[2]  =-10000000.0;

    for (unsigned int i = 0; i < pointCloud->points.size(); i++){

      if ((pointCloud->points[i].x != noDataValue.x) &&
	  (pointCloud->points[i].y != noDataValue.y) &&
	  (pointCloud->points[i].z != noDataValue.z)){
	
	//update point cloud bounding box - START
	if (pointCloud->points[i].x<bbox3.first[0]){
	  bbox3.first[0] = pointCloud->points[i].x;
	}
	if (pointCloud->points[i].x>=bbox3.last[0]){
	  bbox3.last[0]  = pointCloud->points[i].x;
	}
	if (pointCloud->points[i].y<bbox3.first[1]){
	bbox3.first[1] = pointCloud->points[i].y;
	}
	if (pointCloud->points[i].y>=bbox3.last[1]){
	  bbox3.last[1]  = pointCloud->points[i].y;
	}
	if (pointCloud->points[i].z<bbox3.first[2]){
	bbox3.first[2] = pointCloud->points[i].z;
	}
	if (pointCloud->points[i].z>=bbox3.last[2]){
	  bbox3.last[2]  = pointCloud->points[i].z;
	}
      }
    }
    return bbox3;
  }

  //crop a point cloud within a bounding box
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cropToBBOx(pcl::PointCloud<pcl::PointXYZRGB>::Ptr origPC, struct BoundingBox bBox3)
  {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cropPC (new pcl::PointCloud<pcl::PointXYZRGB>);
    cropPC->width=0;
    cropPC->height=1;

    for (unsigned int i = 0; i < origPC->points.size(); i++){ 
      if ((origPC->points[i].x > bBox3.first[0]) && (origPC->points[i].x < bBox3.last[0]) &&
          (origPC->points[i].y > bBox3.first[1]) && (origPC->points[i].y < bBox3.last[1])){
	cropPC->points.push_back(origPC->points[i]);
	cropPC->width = cropPC->width+1;
      }
    }

    cout<<"pc.cc::cropToBBox(): BBOX original PC"<<endl;
    cout<<"pc.cc::cropToBBox(): bBox3.minX="<<bBox3.first[0]<<", bBox3.maxX="<<bBox3.last[0]<<endl;
    cout<<"pc.cc::cropToBBox(): bBox3.minY="<<bBox3.first[1]<<", bBox3.maxY="<<bBox3.last[1]<<endl;
    cout<<"pc.cc::cropToBBox(): bBox3.minZ="<<bBox3.first[2]<<", bBox3.maxZ="<<bBox3.last[2]<<endl;
    
    return cropPC;
  }
#if 0
  //downsamples the point cloud by m_matchPCDownSampleFactor
  //computes the bounding box of the downsampled point cloud
  //removes invalid points (equal to noDataPoint) from the downsampled point cloud
  //returns a new downsampled point cloud
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsamplePC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
						      pcl::PointXYZRGB noDataValue, float downsampleFactor)
  {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr outCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    
    outCloud->width=0;
    outCloud->height=1;

    cout<<"initCloud numPoints = "<<initCloud->points.size()<<", width="<<initCloud->width<<", height="<<initCloud->height<<endl;


    //downsample
    for (unsigned int i=0; i < initCloud->points.size(); i=i+downsampleFactor){
      cout<<"i="<<i<<endl;
      if ((initCloud->points[i].x != noDataValue.x) &&
	  (initCloud->points[i].y != noDataValue.y) &&
	  (initCloud->points[i].z != noDataValue.z) &&
	  (initCloud->points[i].x ==initCloud->points[i].x)&&
	  (initCloud->points[i].y ==initCloud->points[i].y)&&
          (initCloud->points[i].z ==initCloud->points[i].z)) {

	outCloud->points.push_back(initCloud->points[i]);
	outCloud->width = outCloud->width+1;

      }
    }	  

    cout<<"downsamplePC(): numPoints in downsampledPC="<<outCloud->width*outCloud->height<<endl;   
    return outCloud;  
  }
#endif
  //preserves the 2D data structure and does not remove the invalid points
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsamplePC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
						       pcl::PointXYZRGB noDataValue, float downsampleFactor)
  {
    cout<<"pc.cc:: downsamplePC() :initCloud numPoints = "<<initCloud->points.size()
	<<", width="<<initCloud->width<<", height="<<initCloud->height<<endl;

    cout<<"downsampleFactor = "<<downsampleFactor<<endl;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr outCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    
    outCloud->width=0;//initCloud->width/downsampleFactor;
    outCloud->height=1;//initCloud->height/downsampleFactor;

    //downsample
    for (unsigned int row = 0; row < initCloud->height; row = row + downsampleFactor){
      for (unsigned int col = 0; col < initCloud->width; col = col + downsampleFactor){
	int i = row*initCloud->width+col;
	 if ((initCloud->points[i].x != noDataValue.x) &&
	     (initCloud->points[i].y != noDataValue.y) &&
	     (initCloud->points[i].z != noDataValue.z) &&
	     (initCloud->points[i].x ==initCloud->points[i].x)&&
	     (initCloud->points[i].y ==initCloud->points[i].y)&&
             (initCloud->points[i].z ==initCloud->points[i].z)) {

	   outCloud->points.push_back(initCloud->points[row*initCloud->width+col]);
	   outCloud->width = outCloud->width+1;
	   
	 }
      }
    }	  
    
    cout<<"downsamplePC(): numPoints in downsampledPC="<<outCloud->width*outCloud->height<<endl;   
    return outCloud;  
  }
  
  //downsamples the point cloud by m_matchPCDownSampleFactor
  //computes the bounding box of the downsampled point cloud
  //removes invalid points (equal to noDataPoint) from the downsampled point cloud
  //returns a new downsampled point cloud
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr removeNoDataValues(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
						            pcl::PointXYZRGB noDataValue)
  {
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr outCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    outCloud->width=0;
    outCloud->height=1;
 
    //downsample
    for (unsigned int i=0; i < initCloud->points.size(); i++){

      if ((initCloud->points[i].x != noDataValue.x) &&
	  (initCloud->points[i].y != noDataValue.y) &&
	  (initCloud->points[i].z != noDataValue.z) &&
	  (initCloud->points[i].x ==initCloud->points[i].x)&&
	  (initCloud->points[i].y ==initCloud->points[i].y)&&
          (initCloud->points[i].z ==initCloud->points[i].z)) {

	outCloud->points.push_back(initCloud->points[i]);
	outCloud->width = outCloud->width+1;

      }
    }	  

    cout<<"removeNoDataValues(): numPoints in output pointcloud="<<outCloud->width*outCloud->height<<endl;   
    return outCloud;  
  }

  //TODO: move to common/pc.cc
  //this is of more general use than just ICP
  void computeSurfaceNormals(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &points, pcl::PointCloud<pcl::PointNormal>::Ptr &normals, float radiusSearch)
  {
    
     pcl::NormalEstimationOMP<pcl::PointXYZRGB, pcl::PointNormal> n;
     n.setInputCloud(points);

    // Create an empty kdtree representation, and pass it to the normal estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<pcl::PointXYZRGB>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZRGB> ());
    n.setSearchMethod (tree);
    
    // Use all neighbors in a sphere of radius 3cm
    n.setRadiusSearch (radiusSearch);

    pcl::PointCloud<pcl::PointNormal>::Ptr normalsAll(new pcl::PointCloud<pcl::PointNormal>);
  
    //Compute the normals
    n.compute (*normalsAll);

    normals->width = 1;
    normals->height = 0;
    for (int i =0; i < normalsAll->points.size(); i++){
      if (points->points[i].z==points->points[i].z){
        
        pcl::PointNormal thisPointNormal;
	thisPointNormal.x = points->points[i].x; 
	thisPointNormal.y = points->points[i].y; 
	thisPointNormal.z = points->points[i].z;
        thisPointNormal.normal_x = normalsAll->points[i].normal_x;
	thisPointNormal.normal_y = normalsAll->points[i].normal_y;
	thisPointNormal.normal_z = normalsAll->points[i].normal_z;
        normals->points.push_back(thisPointNormal);
        normals->height = (normals->height) + 1;
	
      }
    }
    
    //TODO: remove NAN normals
  }


}
