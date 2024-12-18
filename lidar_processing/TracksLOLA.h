// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef TRACKS_LOLA_H
#define TRACKS_LOLA_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>

//Eigen includes
#include <Eigen/Dense>

//GDAL includes
#include "gdal_priv.h"

using namespace std;

namespace at
{

  struct imgPoint
  {
    double x;
    double y;
    double val;
  };
  
  // regular points on a surface mesh
  struct DEMPoint
  {
    int valid;// map may have holes in the terrain
    double x;
    double y;
    double val;
  };
  
  // from a lidar track, 3D points and time of capture
  class LOLAPoint
  {
  public:
    LOLAPoint();
    LOLAPoint(
	       const Eigen::Vector3f&,//vector<float>, //coords
	       const int          = 0,             // year
	       const int          = 0,             // month
	       const int          = 0,             // day
	       const int          = 0,             // hour
	       const int          = 0,             // min
	       const float        = 0,             // sec
	       const int          = 0              // s
   );
    
  ~LOLAPoint(){};

  double x;
  double y;
  double z;
  int year;
  int month;
  int day;
  int hour; //24hr format
  int min;
  float sec;
  int s; //detector id
  };

  bool isTimeDiff( const LOLAPoint&, const LOLAPoint&, const float);
  
  class LOLAShot {
  public:
    
    explicit LOLAShot( const std::vector<LOLAPoint>& );
    explicit LOLAShot( const LOLAPoint& );
    
    ~LOLAShot(){};
    
    std::vector<LOLAPoint> LOLAPts;
    // points on image and DEM corresponding to point cloud
    std::vector<imgPoint>   imgPt; // original locations of each shot in the original image
    std::vector<DEMPoint>   DEMPt;
    int   valid; //valid shot 1, invalid shot 0. a shot is valid if it has 5 points 
    float reflectance; //the point reflectance for a given camera and light source position; 0 if invalid.
    float image; //the intensity in the original image; -1 if the image is not valid.

    //following variables are for LOLA feature and weight computation
    int   calc_acp;             //is the filter valid here?
    float filter_response;     
    float featurePtRefl;          
    float weightRefl;        
    float featurePtLOLA;
    float filresLOLA;
    

  private:
    void init
      (
       std::vector<LOLAPoint>,
       std::vector<imgPoint> = std::vector<imgPoint>(), // explicitly make empty
       std::vector<DEMPoint> = std::vector<DEMPoint>(), // explicitly make empty
       int              = 0, // valid
       float            = 0, // reflectance
       int              = 0, // calc_acp
       float            = 0, // filter_response
       float            = 0, // featurePtRefl
       float            = 1, // weightRefl
       float            = 0, // featurePtLOLA
       float            = 0  // filresLOLA
       );
  };
  
  LOLAPoint GetPointFromIndex( const LOLAShot& shot, const int );
  int GetPointIndexFromDetector( const LOLAShot& shot, const int detectorIndex );
  
/**
 * Convert a CSV file to a list of LOLA tracks.
 **/
  std::vector<std::vector<LOLAShot> > LOLAFileRead( const std::string& );
  std::vector<std::vector<LOLAShot> > LOLAFilesRead( const std::string& );
  inline std::vector<std::vector<LOLAShot> > CSVFileRead( const std::string& f ){
    // This function is deprecated, please use LOLAFileRead().
    return LOLAFileRead( f );
  }
   
/**
 * Compute reflectance of LOLAShots from camera and sun position, and surface normals.
 **/
  int ComputeAllReflectance( std::vector< std::vector<LOLAShot> >& LOLATracks,
			     const Eigen::Vector3f &cameraPosition,
			     const Eigen::Vector3f &lightPosition);
 
 
  void SaveReflectancePoints( const std::vector<std::vector<LOLAShot> >& LOLATracks,
			      const std::vector<float> &gain_bias,
			      const std::string &filename);
  
  void SaveImagePoints( const std::vector<std::vector<LOLAShot> >& LOLATracks,
			const int detectNum, 
			const std::string& filename);
  
  void SaveAltitudePoints( const std::vector<std::vector<LOLAShot> >& LOLATracks,
			   const int detectNum,
			   const std::string& filename);

  vector<std::string> trackSplitter(std::string lidarFilename, std::string resultsDirname);
  int test(vector<std::string>outTrackFilenames, std::string testTrackListFilename);
  
  void DrawTracksOnImage(vector<vector<LOLAShot> > tracks,
			 GDALDataset* imageDS, string type, int thikness);

  

  
}

#endif /* TRACKS_LOLA_H */

