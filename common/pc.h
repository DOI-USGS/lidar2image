#ifndef PC_H_
#define PC_H_

/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace at
{

  void appendCloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud, std::string pointCloudFilename, int subSampleStep, 
                   int r, int g, int b);
  void appendCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, std::string pointCloudFilename, int subSampleStep);
  void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, std::string filename);
  void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, float xOffset, float yOffset, std::string filename);
  void savePointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, std::string filename);

  void savePointNormal(pcl::PointCloud<pcl::PointNormal> pointCloud, std::string filename);
  void appendPointCloudRGBToFile(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, std::string filename);

  void savePointCloudPtr(pcl::PointCloud<pcl::PointXYZ>::Ptr pointCloud, std::string filename);
  void savePointCloud(pcl::PointCloud<pcl::PointXYZ> pointCloud, std::string filename);
  void saveOrganizedPointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, std::string filename);
  void saveOrganizedPointCloudRGB(pcl::PointCloud<pcl::PointXYZRGB> pointCloud, std::string filename);

  int CountEntries(std::string pointCloudFilename);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr  ReadPCFile(std::string pointCloudFilename, int numPoints, int subSampleStep);
  pcl::PointCloud<pcl::PointXYZ>::Ptr  ReadPCNoRGBFile(std::string pointCloudFilename, int numPoints, int subSampleStep);
  
  //transfor a point cloud using a 4x4 transformation matrix
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr makeTransformedPC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, 
                                                           Eigen::Matrix4f transform);
  struct BoundingBox GetBoundingBox(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pointCloud, pcl::PointXYZRGB noDataValue);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cropToBBOx(pcl::PointCloud<pcl::PointXYZRGB>::Ptr origPC, struct BoundingBox bBox3);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsamplePC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
						      pcl::PointXYZRGB noDataValue, float downsampleFactor);
  //pcl::PointCloud<pcl::PointXYZRGB>::Ptr downsamplePC(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
  //						      float downsampleFactor);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr removeNoDataValues(pcl::PointCloud<pcl::PointXYZRGB>::Ptr initCloud, 
							    pcl::PointXYZRGB noDataValue);
  void computeSurfaceNormals(pcl::PointCloud<pcl::PointXYZRGB>::Ptr &points, pcl::PointCloud<pcl::PointNormal>::Ptr &normals,
			     float radiusSearch);
  
}

#endif	// PC_H_
