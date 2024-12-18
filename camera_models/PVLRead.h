// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

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
