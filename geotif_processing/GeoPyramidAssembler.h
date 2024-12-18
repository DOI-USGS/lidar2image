// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__



#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "gdal_priv.h"

namespace at 
{

    std::vector<struct TerrainFilename> pyrAssembleGeoData(string foreDEMFilename, string backDEMFilename,
							   string foreDRGFilename, string backDRGFilename,
							   struct AssemblerParams assemblerParams, 
							   string resDir, int overwriteFlag);
    
    std::vector<struct TerrainFilename> pyrAssembleGeoData(GDALDataset *foreDEM, GDALDataset *backDEM,
							   GDALDataset *foreDRG, GDALDataset *backDRG,
							   struct bbox unionGeoBBox,
							   float offsetDEM, struct AssemblerParams assemblerParams, 
							   string thisResDir, int overwriteFlag);       
}

