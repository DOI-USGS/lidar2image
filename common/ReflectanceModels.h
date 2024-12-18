// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iomanip>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace at{

  Eigen::Vector3f ComputeNormalFrom3DPointsGeneral(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3); 
  float ComputeLunarLambertianReflectanceFromNormal(Eigen::Vector3f sunPos, Eigen::Vector3f viewPos, 
						    Eigen::Vector3f xyz, Eigen::Vector3f normal); 

}
