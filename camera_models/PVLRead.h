/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#ifndef PVL_READ_H
#define PVL_READ_H

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

//Eigen includes
#include <Eigen/Core>

namespace at
{
  Eigen::Vector3f scanPVLForSunPosition(string pvlFilename);
  Eigen::Vector3f scanPVLForCamPosition(string pvlFilename);
}
#endif//PVL_READ_H
