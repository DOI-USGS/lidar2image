// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../common/StringUtils.h"
#include "../common/Tiling.h"
#include "CoordTransform.h"


using namespace std;

namespace at {

  GDALDatasetH resampleGDALDataset(GDALDataset* origData, float resampleFactor,
				   string resampledFilename, string resizeType, int saveFileFlag);
  GDALDataset* uint16ToF32GDAL(GDALDataset * data, string filename, int saveFileFlag);
  float* uint16ToF32(GDALDataset *data);
  GDALDatasetH resizeGDALDataset(GDALDataset* src, struct bbox geographicBBox, string dstFilename);
  GDALDatasetH cropGDALDataset(GDALDataset* origData,  struct bbox &pixelBBox,
  			       int bandIndex, float upsampleRatio, string filename);
  GDALDatasetH translateGDALDataset(GDALDataset* src, vector<double> newTiePoint,
				    string filename);
  GDALDatasetH changeBandsGDALDataset(GDALDataset* src, int numRasters, string dstFilename);
  void addToGDALDataset(GDALDataset* src, float val);
  float* GetNoDataValue(GDALDataset* data);
  struct bbox Get2DBoundingBox(GDALDataset *data);
  float computeOffset(GDALDataset* foreData, GDALDataset* backData);
  
  std::vector<std::string> findOverlapFiles(std::vector<std::string> tileFilenameArray, GDALDataset* data);
  struct bbox findOverlapPixelBBox(GDALDataset* refDS, GDALDataset *currDS);
  
  vector<double> computeGeographicBBox(GDALDataset* data);
  vector<float>  computeUnionGeographicBBox(vector<string> filenameArray);

  struct bbox changePixelBBoxCoords(GDALDataset* src, GDALDataset* dst, struct bbox srcPixelBBox);  
  struct bbox projectedToPixelBBox(struct bbox projectedBBox, GDALDataset* data);
  struct bbox pixelToProjectedBBox(struct bbox pixelBBox, GDALDataset* data);
  struct bbox geographicToPixelBBox(struct bbox geographicBBox, GDALDataset* data);
  struct bbox pixelToGeographicBBox(struct bbox pixelBBox, GDALDataset* data);
  
  unsigned char* make8UCTile(GDALDataset* ds, float upsampleFactor, 
			     struct TilingParams &tileParams, int bandIndex);
  float*         make32FTile(GDALDataset* ds, float upsampleFactor, 
		             struct TilingParams &tileParams, int bandIndex);
  void   saveGDALDataset(GDALDataset* dataset, float *tileData32F);
  void   saveGDALDataset(GDALDataset* dataset, unsigned char *tileData8UC);
  CPLErr saveGDALDataset(unsigned char *drg, int imageHeight, int imageWidth, int numBands,
			 double*geoTransform, char *pszSRS_WKT_1, std::string filename);
  CPLErr saveGDALDataset(float *dem, int newImageHeight, int newImageWidth,
			 double*geoTransform, char *pszSRS_WKT_1, std::string filename);

}
 
//#endif
