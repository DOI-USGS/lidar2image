// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef PYRTILES_GEOTIF_UTIL_H
#define PYRTILES_GEOTIF_UTIL_H

#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../common/StringUtils.h"
#include "../common/Tiling.h"
#include "GeoUtils.h"
#include "CoordTransform.h"

using namespace std;

typedef struct ResampleParams
{
  int   imageType; //0 (DEM) or 1 (DRG), temporary solution
  int   numPyrLevels; //negative numbers
  float pyrResampleFactor;
  int   imageROILeft;
  int   imageROITop;
  int   imageROIRight;
  int   imageROIBottom;
  int   tileSize;
  int   tilePaddingTop;
  int   tilePaddingLeft;
  int   tilePaddingRight;
  int   tilePaddingBottom;
  int   pyramidMode;
} ResampleParams;

namespace at {

  int  ReadResampleConfigFile(string resampleConfigFilename, struct ResampleParams *resampleParams);
  void PrintResampleParams(struct ResampleParams *resampleParams);
  

  std::vector<std::vector<struct TilingParams> > makeVirtualTiles(GDALDataset *data,
								  /*vector<float> unionGeographicBBox,*/
								  struct bbox unionGeoBBox,
								  struct ResampleParams resampleParams);
  
  std::vector<std::vector<struct TilingParams> > changeVirtualTiles(GDALDataset* demDS, struct ResampleParams demResampleParams,
								    GDALDataset* drgDS, struct ResampleParams drgResampleParams,
								    std::vector<std::vector<struct TilingParams> >pyrTileParamsArray);
  
  void checkVirtualTiles(std::vector<std::vector<struct TilingParams> >pyrTileArray, int tileSize);
 
  void writePyramidTiles(GDALDataset* origDataset, int bandIndex, int numBands,
			 std::vector<std::vector<struct TilingParams> >&pyrTileParamsArray,
			 string mode, string outputDirname);
  
}
 
#endif
