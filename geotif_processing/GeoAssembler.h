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
#include "../common/StringUtils.h"
#include "GeoTiling.h"

namespace at 
{

    struct GeoImage
    {
      GDALDataset* dataset;
      unsigned char* uc8;
      float *f32;
    };

    struct TerrainFilename
    {
      std::vector<string> demFilenames;
      std::vector<string> drgFilenames;
    };

     typedef struct AssemblerParams
    {
      std::vector<float> samplingStep;
      int     weightingMode;
      float   foreMaxPPD;// = 0.03125; 
      int     tileSizeDEM; //128
      std::vector<float> paddingParamsDEM;//0,0,1,1
      int     tileSizeDRG; //512
      std::vector<float> paddingParamsDRG;//1,1,2,2
      float   foreNoDataValDEM;
      float   backNoDataValDEM;
      float   foreNoDataValDRG;
      float   backNoDataValDRG;
      std::vector<float> resampleFactorDRG;
      int     backReference;
      int     pyramidMode;
      std::vector<float> roiDEM;
      std::vector<float> roiDRG; 
    } AssemblerParams;

 
    void saveRegistrationParams(float offset, string filename);


    void PrintSettings(struct AssemblerParams *assemblerParams);
    int ReadAssemblerConfigFile(std::string assemblerConfigFilename, struct AssemblerParams *assemblerParams);
    
    GDALDatasetH makeReferenceLayer(GDALDataset *fore, GDALDataset *back,
				  struct bbox unionGeoBBox,
				  struct AssemblerParams assemblerParams,
				  string mode, float*foreToReferenceResampleFactor);
  
    vector<vector<string> > assembleGeoData(string foreFilename, string backFilename,
				   struct AssemblerParams assemblerParams, 
                                   string mode, string resDir, int overwriteFlag);
    
    vector<vector<string> > assembleGeoData(GDALDataset* fore,  GDALDataset* back,
					    struct AssemblerParams assemblerParams, 
					    std::string mode,
					    struct bbox unionGeoBBox,
					    std::string foreFilename, std::string resDir, int overwriteFlag);

    
    vector<struct TerrainFilename> assembleGeoData(string foreDEMFilename, string backDEMFilename,
						   string foreDRGFilename, string backDRGFilename,
						   struct AssemblerParams assemblerParams, 
						   string resDir, int overwriteFlag);
    
    vector<struct TerrainFilename>assembleGeoData(GDALDataset *foreDEM, GDALDataset *backDEM,
						  GDALDataset *foreDRG, GDALDataset *backDRG,
						  vector<float> unionBBoxGeographic,
						  float offsetDEM, struct AssemblerParams assemblerParams, 
						  string thisResDir, int overwriteFlag);
     
    std::vector<std::vector<std::string> > writeAssembledTiles(GDALDataset* fore, GDALDataset* back, GDALDatasetH hReferenceDS,
							       float offsetDEM, int weightingMode,float backUpsamplingFactor,
							       int overwriteFlag, int reference,
							       std::vector<std::vector<struct TilingParams> > pyrTileParamsArray,
							       string mode, string thisResDir);
      
    vector<float> computeResampleFactors(GDALDataset* foreData, GDALDataset* backData,
					  int backReference, float foreMaxMPP);
   
}

