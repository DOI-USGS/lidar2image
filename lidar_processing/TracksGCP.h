// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef GCP_H
#define GCP_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>

//Eigen includes
#include <Eigen/Dense>

namespace at
{

  struct gcp{
    float lon;
    float lat;
    float rad;
    int   year;
    int   month;
    int   day;
    int   hour;
    int   min;
    float sec;

    float sigma_lon;
    float sigma_lat;
    float sigma_rad;
    std::vector<std::string> filename;
    int trackIndex;
    int shotIndex;
    vector<float> x;
    vector<float> y;
    vector<float> x_before;
    vector<float> y_before;
  };

  vector<gcp> ComputeSalientLOLAFeatures(std::vector<std::vector<LOLAShot> > & trackPts);

  void SaveGCPoints( const std::vector<gcp>&, const std::string& );
  void SaveGCP( const vector<gcp>& gcpArray, const string& gcpFilename);

}
#endif /* GCP_H */

