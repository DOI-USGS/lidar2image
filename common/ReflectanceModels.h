/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <iomanip>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace at{

  Eigen::Vector3f ComputeNormalFrom3DPointsGeneral(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3); 
  float ComputeLunarLambertianReflectanceFromNormal(Eigen::Vector3f sunPos, Eigen::Vector3f viewPos, 
						    Eigen::Vector3f xyz, Eigen::Vector3f normal); 

}
