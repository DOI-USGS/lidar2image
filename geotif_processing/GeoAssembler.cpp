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
#include "gdalwarper.h"
#include <ogr_spatialref.h>

#include "CoordTransform.h"
#include "GeoAssembler.h"
#include "GeoUtils.h"
#include "GeoTiling.h"

#include "../common/StringUtils.h"
#include "../common/FileListUtils.h"
#include "../common/ImageProcessing.h"
#include "../common/FileSystemUtils.h"

using namespace std;

namespace at 
{

  inline bool IsFloatNoData(float value, double noDataValue)
  {
    if (isnan(noDataValue))
      return isnan(value);

    return (noDataValue == value);
  }
    
  void saveRegistrationParams(float offset, string filename)
  {
    ofstream outfile(filename.c_str());
    outfile<<"ROTATION_MATRIX: 1 0 0 0 1 0 0 0 1"<<endl;
    outfile<<"TRANSLATION_VECTOR: 0 0 0"<<endl;
    outfile<<"LON_LAT_RAD_OFFSET: 0 0 "<<offset<<endl;
    outfile.close();
  }
  
  void PrintSettings(struct AssemblerParams *assemblerParams)
  {
    
    cout<<"WEIGHTING_MODE: "<<assemblerParams->weightingMode<<endl;
    
    cout<<"FORE_MAX_PPD: "<<assemblerParams->foreMaxPPD<<endl;
    cout<<"SAMPLING_STEP: "<<assemblerParams->samplingStep[0]<<", "<<assemblerParams->samplingStep[0]<<endl;
    
    cout<<"TILE_SIZE_DEM: "<<assemblerParams->tileSizeDEM<<endl;
    cout<<"PADDING_PARAMS_DEM: "<<assemblerParams->paddingParamsDEM[0]<<", "<<assemblerParams->paddingParamsDEM[1]<<", "
	<<assemblerParams->paddingParamsDEM[2]<<", "<<assemblerParams->paddingParamsDEM[3]<<endl;
    cout<<"TILE_SIZE_DRG: "<<assemblerParams->tileSizeDRG<<endl;
    cout<<"PADDING_PARAMS_DRG: "<<assemblerParams->paddingParamsDRG[0]<<", "<<assemblerParams->paddingParamsDRG[1]<<", "
	<<assemblerParams->paddingParamsDRG[2]<<", "<<assemblerParams->paddingParamsDRG[3]<<", "<<endl;
    cout<<"FORE_NO_DATA_VAL_DEM: "<<assemblerParams->foreNoDataValDEM<<endl;
    cout<<"BACK_NO_DATA_VAL_DEM: "<<assemblerParams->backNoDataValDEM<<endl;
    cout<<"FORE_NO_DATA_VAL_DRG :"<<assemblerParams->foreNoDataValDRG<<endl;
    cout<<"BACK_NO_DATA_VAL_DRG: "<<assemblerParams->backNoDataValDRG<<endl;
    cout<<"BACK_REFERENCE: "<<assemblerParams->backReference<<endl;
    cout<<"PYRAMID_MODE: "<<assemblerParams->pyramidMode<<endl;
    cout<<"ROI_DEM: "<<assemblerParams->roiDEM[0]<<", "<<assemblerParams->roiDEM[1]<<", "
	<<assemblerParams->roiDEM[2]<<", "<<assemblerParams->roiDEM[3]<<endl;
    cout<<"ROI_DRG: "<<assemblerParams->roiDRG[0]<<", "<<assemblerParams->roiDRG[1]<<", "
	<<assemblerParams->roiDRG[2]<<", "<<assemblerParams->roiDRG[3]<<", "<<endl;
    
    cout<<"RESAMPLE_FACTOR_DRG: "<<assemblerParams->resampleFactorDRG[0]<<", "<<assemblerParams->resampleFactorDRG[1]<<endl;
  }
  
  int ReadAssemblerConfigFile(string assemblerConfigFilename, struct AssemblerParams *assemblerParams)
  {
    ifstream configFile (assemblerConfigFilename.c_str());
    std::string line;
    std::string identifier;
    
    //matching params - START
    assemblerParams->weightingMode = 0;
    
    assemblerParams->samplingStep.resize(2);
    assemblerParams->samplingStep[0] = 256;
    assemblerParams->samplingStep[1] = 256;
    
    //tiling params - START
    assemblerParams->tileSizeDEM=128;
    assemblerParams->paddingParamsDEM.resize(4);
    assemblerParams->paddingParamsDEM[0]=0;
    assemblerParams->paddingParamsDEM[1]=0;
    assemblerParams->paddingParamsDEM[2]=1;
    assemblerParams->paddingParamsDEM[3]=1;
    
    assemblerParams->tileSizeDRG=512;
    assemblerParams->paddingParamsDRG.resize(4);
    assemblerParams->paddingParamsDRG[0]=1;
    assemblerParams->paddingParamsDRG[1]=1;
    assemblerParams->paddingParamsDRG[2]=2;
    assemblerParams->paddingParamsDRG[3]=2;
    
    assemblerParams->foreNoDataValDEM = -3.4028226550889e+38;
    assemblerParams->backNoDataValDEM = -3.4028226550889e+38;
    assemblerParams->foreNoDataValDRG = 0;
    assemblerParams->backNoDataValDRG = 0;
    
    assemblerParams->backReference = 1;
    assemblerParams->pyramidMode = 1;
    
    assemblerParams->roiDEM.resize(4);
    assemblerParams->roiDEM[0]=0;
    assemblerParams->roiDEM[1]=0;
    assemblerParams->roiDEM[2]=0;
    assemblerParams->roiDEM[3]=0;
    
    assemblerParams->roiDRG.resize(4);
    assemblerParams->roiDRG[0]=0;
    assemblerParams->roiDRG[1]=0;
    assemblerParams->roiDRG[2]=0;
    assemblerParams->roiDRG[3]=0;
    
    assemblerParams->resampleFactorDRG.resize(2);
    assemblerParams->resampleFactorDRG[0]=4.0;
    assemblerParams->resampleFactorDRG[1]=4.0;
    //matching params - END
    
    if (configFile.is_open()){ 
      
      cout << "Reading settings file "<< assemblerConfigFilename<<"."<<endl; 
      
      while (!configFile.eof()){
	std::getline(configFile, line);
	
	if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines
	  
	  stringstream lineStream(line);
	  
	  string keyword("");
	  lineStream >> keyword;
	  
	  if (keyword.compare(string("WEIGHTING_MODE"))==0){
	    string stringVal; 
	    lineStream >> stringVal;
	    assemblerParams->weightingMode = (int)floor(atof(stringVal.c_str()));
	  }
	  
	  if (keyword.compare(string("FORE_MAX_PPD"))==0){
	    string val1; 
	    lineStream >> val1;
	    assemblerParams->foreMaxPPD = atof(val1.c_str());
	  }
	  if (keyword.compare(string("SAMPLING_STEP"))==0){
	    string val1, val2; 
	    lineStream >> val1>>val2;
	    assemblerParams->samplingStep[0]=(int)floor(atof(val1.c_str()));
	    assemblerParams->samplingStep[1]=(int)floor(atof(val2.c_str()));
	  }
	  
	  if (keyword.compare(string("TILE_SIZE_DEM"))==0){
	    string val1; 
	    lineStream >> val1;
	    assemblerParams->tileSizeDEM=(int)floor(atof(val1.c_str()));
	  }
	  if (keyword.compare(string("PADDING_PARAMS_DEM"))==0){
	    string val1, val2, val3, val4; 
	    lineStream >> val1>>val2>>val3>>val4;
	    assemblerParams->paddingParamsDEM[0]=(int)floor(atof(val1.c_str()));
	    assemblerParams->paddingParamsDEM[1]=(int)floor(atof(val2.c_str()));
	    assemblerParams->paddingParamsDEM[2]=(int)floor(atof(val3.c_str()));
	    assemblerParams->paddingParamsDEM[3]=(int)floor(atof(val4.c_str()));
	  }
	  if (keyword.compare(string("TILE_SIZE_DRG"))==0){
	    int val1; 
	    lineStream >> val1;
	    assemblerParams->tileSizeDRG=val1;
	  }
	  if (keyword.compare(string("PADDING_PARAMS_DRG"))==0){
	    int val1, val2, val3, val4; 
	    lineStream >> val1>>val2>>val3>>val4;
	    assemblerParams->paddingParamsDRG[0]=val1;
	    assemblerParams->paddingParamsDRG[1]=val2;
	    assemblerParams->paddingParamsDRG[2]=val3;
	    assemblerParams->paddingParamsDRG[3]=val4;
	  }
	  
	  if (keyword.compare(string("FORE_NO_DATA_VAL_DEM"))==0){
	    float val1; 
	    lineStream >> val1;
	    assemblerParams->foreNoDataValDEM=val1;
	  }
	  
	  if (keyword.compare(string("BACK_NO_DATA_VAL_DEM"))==0){
	    float val1; 
	    lineStream >> val1;
	    assemblerParams->backNoDataValDEM=val1;
	  }
	  
	  if (keyword.compare(string("FORE_NO_DATA_VAL_DRG"))==0){
	    float val1; 
	    lineStream >> val1;
	    assemblerParams->foreNoDataValDRG=val1;
	  }
	  
	  if (keyword.compare(string("BACK_NO_DATA_VAL_DRG"))==0){
	    float val1; 
	    lineStream >> val1;
	    assemblerParams->backNoDataValDRG=val1;
	  }
	  
	  if (keyword.compare(string("BACK_REFERENCE"))==0){
	    int val1; 
	    lineStream >> val1;
	    assemblerParams->backReference=val1;
	  }
	  
	  if (keyword.compare(string("PYRAMID_MODE"))==0){
	    int val1; 
	    lineStream >> val1;
	    assemblerParams->pyramidMode=val1;
	  }
	  
	  if (keyword.compare(string("ROI_DEM"))==0){
	    int val1, val2, val3, val4; 
	    lineStream >> val1>>val2>>val3>>val4;
	    assemblerParams->roiDEM[0]=val1;
	    assemblerParams->roiDEM[1]=val2;
	    assemblerParams->roiDEM[2]=val3;
	    assemblerParams->roiDEM[3]=val4;
	  }
	  if (keyword.compare(string("ROI_DRG"))==0){
	    int val1, val2, val3, val4; 
	    lineStream >> val1>>val2>>val3>>val4;
	    assemblerParams->roiDRG[0]=val1;
	    assemblerParams->roiDRG[1]=val2;
	    assemblerParams->roiDRG[2]=val3;
	    assemblerParams->roiDRG[3]=val4;
	  }
	  
	  if (keyword.compare(string("RESAMPLE_FACTOR_DRG"))==0){
	    float val1, val2; 
	    lineStream >> val1>>val2;
	    assemblerParams->resampleFactorDRG[0]=val1;
	    assemblerParams->resampleFactorDRG[1]=val2;
	  }
	  
	}
	
      }
      
      configFile.close();
    }
    
    return 0;
  }
  
  GDALDatasetH makeReferenceLayer(GDALDataset *fore, GDALDataset *back,
				  struct bbox unionGeoBBox,
				  struct AssemblerParams assemblerParams,
				  string mode, float*foreToReferenceResampleFactor)
  {

    string refScaledFilename;
    string refFilename;
    
    if (mode.compare("DEM")==0){
      refScaledFilename="refDEMScaled.tiff";
      refFilename="refDEM.tiff";
    }
    if (mode.compare("DRG")==0){
      refScaledFilename="refDRGScaled.tiff";
      refFilename="refDRG.tiff";
    }
    vector<float> upsampleFactors;
    float minMPP;
    minMPP=assemblerParams.foreMaxPPD;
    if (mode.compare("DRG")==0){
      minMPP =assemblerParams.foreMaxPPD/4.0;
    }
    upsampleFactors = computeResampleFactors(fore, back,
					     assemblerParams.backReference,
					     minMPP);

    GDALDatasetH refScaledH;
    int numForeBands = fore->GetRasterCount();
    int numBackBands = back->GetRasterCount();
    int numRefBands=numBackBands;
    if (numForeBands>numRefBands){numRefBands=numForeBands;}
    cout<<"assembler_util.cc:: makeReferenceLayer(): numForeBands="<<numForeBands
	<<", numBackBands="<<numBackBands<<", numRefBands="<<numRefBands<<endl;

    GDALDatasetH refLayerH ;
    
    if ((assemblerParams.backReference==1)||(assemblerParams.backReference==3)){//back is reference
      //reference uses the sampling grid of background and background resolution
      cout<<"makeReferenceLayer(): BACK reference: back upsampleFactor"<< upsampleFactors[0]<<", fore upsampleFactor "<<upsampleFactors[1]<<endl;
      
      refScaledH = resampleGDALDataset(back, 1.0/upsampleFactors[0], refScaledFilename, string("round"), 0);
      *foreToReferenceResampleFactor = 1.0/upsampleFactors[1];

      GDALDatasetH refScaledMultiBandH = changeBandsGDALDataset((GDALDataset*)refScaledH, numRefBands, "multiband.tiff");
      int numRefTempBands = ((GDALDataset*)refScaledMultiBandH)->GetRasterCount();
      
      bbox refDEMPixelBBox = geographicToPixelBBox(unionGeoBBox, (GDALDataset*)refScaledH);
       
      //generate a reference DEM (reference DEM resolution, unionBBox coverage) 
      refLayerH = resizeGDALDataset((GDALDataset*)refScaledMultiBandH, refDEMPixelBBox, refFilename);

      GDALClose(refScaledH);
      GDALClose(refScaledMultiBandH);
    }
    else{//fore is reference
      //reference DEM uses the sampling grid of foreground and background resolution
      cout<<"makeReferenceLayer(): FORE reference: back upsampleFactor: "<< upsampleFactors[0]<<", fore upsampleFactor: "<<upsampleFactors[1]<<endl;
      
      refScaledH = resampleGDALDataset(fore, 1.0/upsampleFactors[0], refScaledFilename, string("round"),  0);
      
      //and after resampling
      *foreToReferenceResampleFactor = upsampleFactors[0];

      GDALDatasetH refScaledMultiBandH = changeBandsGDALDataset((GDALDataset*)refScaledH, numRefBands, "multiband.tiff");
      int numRefTempBands = ((GDALDataset*)refScaledMultiBandH)->GetRasterCount();
            
      bbox refDEMPixelBBox = geographicToPixelBBox(unionGeoBBox, (GDALDataset*)refScaledH);

      refDEMPixelBBox.xl = round(refDEMPixelBBox.xl/(int)upsampleFactors[0])*upsampleFactors[0];
      refDEMPixelBBox.yt = round(refDEMPixelBBox.yt/(int)upsampleFactors[0])*upsampleFactors[0];
      refDEMPixelBBox.xr = round(refDEMPixelBBox.xr/(int)upsampleFactors[0])*upsampleFactors[0];
      refDEMPixelBBox.yb = round(refDEMPixelBBox.yb/(int)upsampleFactors[0])*upsampleFactors[0];
      /*
      refDEMPixelBBox.xl = round(refDEMPixelBBox.xl/16)*16;
      refDEMPixelBBox.yt = round(refDEMPixelBBox.yt/16)*16;
      refDEMPixelBBox.xr = round(refDEMPixelBBox.xr/16)*16;
      refDEMPixelBBox.yb = round(refDEMPixelBBox.yb/16)*16;
      */
      //generate a reference DEM (reference DEM resolution, unionBBox coverage) 
      refLayerH = resizeGDALDataset((GDALDataset*)refScaledMultiBandH, refDEMPixelBBox, refFilename);

      GDALClose(refScaledH);
      GDALClose(refScaledMultiBandH);
    
    }

  #if 0
    double adfGeoTransform[6];
    ((GDALDataset*)refLayerH)->GetGeoTransform(adfGeoTransform);
    cout<<"assembler_util.cc:: makeReferenceLayer(): refLayer pixel resolution="
	<<adfGeoTransform[1]<<", "<<adfGeoTransform[5]<<endl;
    vector<double> projected;
    projected.resize(2);
    //HIRISE
    projected[0]= 8131125.0;
    projected[1]=-259188.0;
    vector<double> pixel = projectedToPixel(projected, (GDALDataset*)refLayerH/*refScaledMultiBandH*/);
    cout<<"assembler_util.cc:: makeReferenceLayer(): HIRISE tiepoint pixel="<<pixel[0]<<", "<<pixel[1]<<endl;
    
    //HRSC
    projected[0]= 8066131.994;
    projected[1]=-200891.239;
    //vector<double>
    pixel = projectedToPixel(projected, (GDALDataset*)refLayerH/*refScaledMultiBandH*/);
    cout<<"assembler_util.cc:: makeReferenceLayer(): HRSC tiepoint pixel="<<pixel[0]<<", "<<pixel[1]<<endl;
    
    //vector<double> pixel;
    pixel.resize(2);    
    pixel[0]=0; pixel[1]=0;
    //vector<double>
    projected = pixelToProjected(pixel, (GDALDataset*)refLayerH);
    cout<<"assembler_util.cc:: makeReferenceLayer(): BackTiePoint="<<projected[0]<<", "<<projected[1]<<endl;
 #endif


    //exit(1);
    
    return refLayerH;
    
  }
  
  //compute fore, back resampling factors
  vector<float> computeResampleFactors(GDALDataset* foreData1, GDALDataset* backData1,
				       int backReference, float foreMaxMPP)
  {
    double backGeoTransform[6];
    backData1->GetGeoTransform(backGeoTransform);
    
    //remove the ROI from the back image, added by Ara , Dec 6th, 2014
    //backGeoTransform[0]=backGeoTransform[0]-roi[0]*backGeoTransform[1];
    //backGeoTransform[3]=backGeoTransform[3]-roi[1]*backGeoTransform[5];
    
    double foreGeoTransform[6];
    foreData1->GetGeoTransform(foreGeoTransform);
    float initForeOverBackResolutionRatio = fabs(backGeoTransform[1]/foreGeoTransform[1]);
    
    //adjust the upsampling factors for foreground and background in case we use a quadtree pyramid
    
    int numSubPyrLevels;
    float backUpsamplingFactor;
    float foreUpsamplingFactor;
    
    //return numSubPyrLevels, backUpsamplingFactor, foreUpsamplingFactor 
    cout<<"assembler_util.cc:: computeResampleFactors(): backReference="<<backReference<<endl;
    switch(backReference)
      {
      case 0: //foreground is the reference quadtree
	{        
	  cout<<"assembler_util.cc:: computeResampleFactors(): FORE Reference in quadtree"<<endl;
	  foreUpsamplingFactor = 1;
	  //compute the backUpsamplingFactor and the number of levels corresponding to a power-of-two pyramid
	  numSubPyrLevels = (int)ceil(log2(initForeOverBackResolutionRatio));
	  cout<<"assembler_util.cc:: computeResampleFactors(): numSubPyrLevels="<<numSubPyrLevels<<endl;
	  backUpsamplingFactor = pow(2, (float)numSubPyrLevels);
	  cout<<"assembler_util.cc:: computeResampleFactors(): backUpsamplingFactor="<<backUpsamplingFactor<<endl;
	}
	break;
	
      case 1: //background is reference in quadtree
	{
	  cout<<"assembler_util.cc:: computeResampleFactors(): BACK Reference in quadtree"<<endl;
	  backUpsamplingFactor = 1;
	  //set a max value for the foreUpsamplingFactor;
	  cout<<"assembler_util.cc:: computeResampleFactors(): foreGeoTransform[1]="<<foreGeoTransform[1]<<endl;
	  cout<<"assembler_util.cc:: computeResampleFactors(): foreMaxMPP="<<foreMaxMPP<<endl;
	  if (foreGeoTransform[1] < foreMaxMPP){
	    foreGeoTransform[1]=foreMaxMPP;
	  }
	  foreUpsamplingFactor = foreGeoTransform[1]/backGeoTransform[1];
	  //compute the initial resolution ratio between foreground and background
	  float foreOverBackResolutionRatio = fabs(backGeoTransform[1]/foreGeoTransform[1]);
	  numSubPyrLevels = (int)ceil(log2(foreOverBackResolutionRatio));
          
	  cout<<"assembler_util.c:: computeResampleFactors():: foreOverBackResolutionRatio Init="<<foreOverBackResolutionRatio<<endl;
	}
	break;
      case 2:// foreground is reference, no quadtree
	{
	  cout<<"assembler_util.cc:: computeResampleFactors: FORE Reference"<<endl;
	  numSubPyrLevels = 1;
	  foreUpsamplingFactor = 1;
	  backUpsamplingFactor = initForeOverBackResolutionRatio;
	}
	break;
      case 3: //background is reference, no quadtree
	{
	  cout<<"assembler_util.cc:: computeResampleFactors(): BACK Reference"<<endl;
	  numSubPyrLevels = 1;
	  backUpsamplingFactor = 1;
	  foreUpsamplingFactor = 1/initForeOverBackResolutionRatio;
	}
	break;
      default://default do nothing: should print invalid option and exit with an error.
	break;
      }  
    
    cout<<"assembler_util.cc:: computeResampleFactors(): numSubPyrLevels="<<numSubPyrLevels<<endl;
    cout<<"assembler_util.cc:: computeResampleFactors(): foreUpsamplingFactor="<<foreUpsamplingFactor
	<<", backUpsamplingFactor="<<backUpsamplingFactor<<endl;   
    cout<<"---------------------------------"<<endl;
    
    vector<float> upsampleFactors; upsampleFactors.resize(2);
    upsampleFactors[0]= backUpsamplingFactor;
    upsampleFactors[1]= foreUpsamplingFactor;
    return upsampleFactors;
  }
  
  //tileParams is given in the original size background image
  //upsampleRatioBackImag is the upsampling factor of the original background
  //image to match the asembled image
  //upsampleRatioForeImg is the upsampling factor of the original foreground
  //image to match the asembled image
  //offset is the elevation offset
   struct GeoImage writeAssembledTile(GDALDataset* foreDS, GDALDataset* backDS, GDALDatasetH hReferenceDS,
				      float elevationOffset, int weightingMode,
				      float foreToReferenceResampleFactor, 
				      struct bbox &tilePixelBBoxRef, string fileType,
				      int bandIndex, int backReference, string assembledImgFilename)
  {
        
    struct GeoImage geoImage;
    geoImage.uc8 = NULL;
    geoImage.f32 = NULL;
    geoImage.dataset = NULL;
	    
    //determine the tile filenames - START
    string  assembledAccFilename; 
    string  assembledPCFilename;
    string  foregroundImgFilename;
 
    if (fileType.compare("DEM")==0){
      assembledAccFilename = assembledImgFilename;
      FindAndReplace(assembledAccFilename, "_dem", "_acc");
      
      assembledPCFilename = assembledImgFilename;
      FindAndReplace(assembledPCFilename, "_dem.tif", "_pc.txt");
      
      foregroundImgFilename = assembledImgFilename;
      FindAndReplace(foregroundImgFilename, "assembled_", "fore_");
    }
    if (fileType.compare("DRG")==0){
      foregroundImgFilename = assembledImgFilename;
      FindAndReplace(foregroundImgFilename, "assembled_", "fore_");
    }

    //determine the tile filenames - END
    struct bbox backBBox;
    backBBox.xl=0;
    backBBox.yt=0;
    backBBox.xr=GDALGetRasterXSize(backDS);
    backBBox.yb=GDALGetRasterYSize(backDS);
  
    struct bbox foreBBox;
    foreBBox.xl=0;
    foreBBox.yt=0;
    foreBBox.xr=GDALGetRasterXSize(foreDS);
    foreBBox.yb=GDALGetRasterYSize(foreDS);
 
    struct bbox tileGeoBBox;    
    tileGeoBBox = pixelToGeographicBBox(tilePixelBBoxRef, (GDALDataset*)hReferenceDS);

    //compute pixel coordinates of the back/tile intersection in the back - START
    struct bbox backTilePixelBBoxBackTemp = geographicToPixelBBox(tileGeoBBox, backDS);
    struct bbox backTilePixelBBoxBack = computeBoxIntersection(backTilePixelBBoxBackTemp, backBBox);
    //compute pixel coordinates of the back/tile intersection in the back image - START

    //compute pixel coordinates of the fore/tile intersection in the fore - START
    struct bbox foreTilePixelBBoxForeTemp = geographicToPixelBBox(tileGeoBBox, foreDS);
    struct bbox foreTilePixelBBoxFore = computeBoxIntersection(foreTilePixelBBoxForeTemp, foreBBox);
    //compute pixel coordinates of the tile in the fore image - END

    //compute pixel coordinates of the fore/tile intersection in the (back)reference - START
    struct bbox foreTilePixelBBoxRefTemp = changePixelBBoxCoords(foreDS, (GDALDataset*)hReferenceDS, foreTilePixelBBoxFore);
    struct bbox foreTilePixelBBoxRef = computeBoxIntersection(foreTilePixelBBoxRefTemp, tilePixelBBoxRef);
    //compute pixel coordinates of the fore/tile intersection in the fore image - END
    
    //compute pixel coordinates of the back/tile intersection in the reference - START
    struct bbox backTilePixelBBoxRefTemp = changePixelBBoxCoords(backDS, (GDALDataset*)hReferenceDS, backTilePixelBBoxBack);
    struct bbox backTilePixelBBoxRef = computeBoxIntersection(backTilePixelBBoxRefTemp, tilePixelBBoxRef);
    //compute pixel coordinates of the fore/tile intersection in the fore image - END

    if ((isRectangle(foreTilePixelBBoxFore)) && (isRectangle(backTilePixelBBoxBack))){
      
      //create a back and tile GDALDataset
      GDALDataset* backTileDS;
      GDALDatasetH hBackTileDS;
      
      hBackTileDS = GDALOpen(assembledImgFilename.c_str(), GA_Update); 
      
      if (((GDALDataset*)hBackTileDS)==NULL){
	GDALDatasetH hBackTileTempDS = cropGDALDataset((GDALDataset*)hReferenceDS, tilePixelBBoxRef, bandIndex,
						       foreToReferenceResampleFactor, "tile.tif");
	hBackTileDS = changeBandsGDALDataset((GDALDataset*)hBackTileTempDS, 1, assembledImgFilename);
	GDALClose (hBackTileTempDS);
	
	//reproject the backData on the tile
	CPLErr err = GDALReprojectImage (backDS, backDS->GetProjectionRef(), hBackTileDS, NULL,
					 GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL); 		
      }
      backTileDS = (GDALDataset*) hBackTileDS;
      //create a tile at the fore resolution - END
      
      int backTileWidthBack = GDALGetRasterXSize(backTileDS);
      int backTileHeightBack = GDALGetRasterYSize(backTileDS);
      cout<<"assembler_util.cc:: writeTile(): backTileWidthBack="<<backTileWidthBack
	  <<", backTileHeightBack="<<backTileHeightBack<<endl;
      
      //access the projected background image on the tile - START 
      if ((backTileWidthBack>0)&&(backTileHeightBack>0)){
	
	float* tileData32F = NULL;
	unsigned char *tileData8UC = NULL;

	GDALRasterBand* tileBand = backTileDS->GetRasterBand(1);
	GDALDataType tileDataType = tileBand->GetRasterDataType();
      
	if (tileDataType==GDT_Byte){	
           tileData8UC = new unsigned char[backTileWidthBack*backTileHeightBack];
	   cout<<"assembler_util.cc:: writeAssembledTile(): Back is BYTE"<<endl;
	   tileBand->RasterIO(GF_Read, 0, 0, backTileWidthBack, backTileHeightBack,
	                      tileData8UC, backTileWidthBack, backTileHeightBack, 
	                      tileDataType, 0, 0);
	   
        }
	if (tileDataType==GDT_Float32){
            tileData32F = new float[backTileWidthBack*backTileHeightBack];
	    cout<<"assembler_util.cc:: writeAssembledTile(): Back is FLOAT"<<endl;
	    tileBand->RasterIO(GF_Read, 0, 0, backTileWidthBack, backTileHeightBack,
	                       tileData32F, backTileWidthBack, backTileHeightBack, 
			       tileDataType, 0, 0);	
        }
        //access the projected background image on the tile - END
     
     
        //create a fore and tile GDALDataset
	GDALDatasetH hForeTileDS;
	
	hForeTileDS = cropGDALDataset((GDALDataset*)hReferenceDS, foreTilePixelBBoxRef, bandIndex,
				      foreToReferenceResampleFactor, foregroundImgFilename);
	
	//reproject the foreground on the fore and tile GDALdataset
	CPLErr err_1 = GDALReprojectImage (foreDS, foreDS->GetProjectionRef(), hForeTileDS, NULL,
					   GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL); 
	GDALDataset *foreTileDS = (GDALDataset*) hForeTileDS;
	
	int foreTileWidthBack;
	int foreTileHeightBack;
	
	foreTileWidthBack = GDALGetRasterXSize(foreTileDS);
	foreTileHeightBack = GDALGetRasterYSize(foreTileDS);
		  
	if ((foreTileWidthBack>0) && (foreTileHeightBack>0)){

	  float * foreTileData32F= NULL; //reprojected fore data 32F over this background tile
	  unsigned char *foreTileData8UC = NULL; //reprojected fore data 8UC over this background tile
	  
	  //access the projected foreground image on the tile - START

	  int thisNumBands=foreTileDS->GetRasterCount();
	  GDALRasterBand* foreTileBand = foreTileDS->GetRasterBand(bandIndex);
	  
	  GDALDataType foreTileDataType=foreTileBand->GetRasterDataType();
	  
	  if (foreTileDataType==GDT_Byte){ 
	    foreTileData8UC = new unsigned char[foreTileHeightBack*foreTileWidthBack];
	    cout<<"assembler_util.cc:: writeAssembledTile(): Fore is BYTE"<<endl;
	    //read foreBand and copy to rForeData8UC
	    foreTileBand->RasterIO(GF_Read, 0,0,
				   foreTileWidthBack, foreTileHeightBack,
			           foreTileData8UC, foreTileWidthBack, foreTileHeightBack, 
				   foreTileDataType, 0, 0);
	    cout<<"assembler_util.cc:: writeAssembledTile(): Fore is BYTE after"<<endl;
	  }
	  if (foreTileDataType==GDT_Float32){
	    foreTileData32F = new float[foreTileHeightBack*foreTileWidthBack];
	    cout<<"assembler_util.cc:: writeAssembledTile(): Fore is FLOAT"<<endl;
	    //read foreBand and copy to rForeData32F
	    foreTileBand->RasterIO(GF_Read, 0,0,
				   foreTileWidthBack, foreTileHeightBack, 
				   foreTileData32F, foreTileWidthBack, foreTileHeightBack, 
				   foreTileDataType, 0, 0);
	
	    holeFilling(foreTileData32F, foreTileWidthBack, foreTileHeightBack, -100000.0);
	  }
	  //access the projected foreground image on the tile - END
	  
	  //fill in the foreground values - START      
	  double foreNoDataVal, foreDSNoDataVal;
	  int pbSuccess = 0;
	  
	  if (foreTileDataType == GDT_Byte){
	    foreNoDataVal = 0.0;
	    foreDSNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)foreTileBand, &pbSuccess);
	    if (pbSuccess == 1) {
	      foreNoDataVal = foreDSNoDataVal;
	    } 
	  }
	  else if (foreTileDataType == GDT_Float32) {
	    foreNoDataVal = NAN;
	    foreDSNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)foreTileBand, &pbSuccess);
	    if (pbSuccess == 1) {
	      foreNoDataVal = foreDSNoDataVal;
	    } 
	  }
	  cout<<"assembler_util.cc:: writeAssembledTile(): foreNoDataVal="<<foreNoDataVal<<endl;
	  
	  vector<double> newForePixel;
	  newForePixel.resize(2);
	  newForePixel[0]=0;
	  newForePixel[1]=0;     
	  
	  int foreStartRow = (round(foreTilePixelBBoxRef.yt)-tilePixelBBoxRef.yt)*foreToReferenceResampleFactor;
	  int foreEndRow   = (round(foreTilePixelBBoxRef.yb)-tilePixelBBoxRef.yt)*foreToReferenceResampleFactor;
	  int foreStartCol = (round(foreTilePixelBBoxRef.xl)-tilePixelBBoxRef.xl)*foreToReferenceResampleFactor;
	  int foreEndCol   = (round(foreTilePixelBBoxRef.xr)-tilePixelBBoxRef.xl)*foreToReferenceResampleFactor;
	  	  
	  //for every pixel in the foreground region of the assembled (and upsampled) tile
	  for (int j = foreStartRow; j < foreEndRow; j++){
	    newForePixel[0]=0;
	    for (int i = foreStartCol; i < foreEndCol; i++){	  
	      //check not to go outside the boundaries of the foreground image
	      if ((newForePixel[0] >=0) && (newForePixel[0] < foreTileWidthBack) &&
		  (newForePixel[1] >=0) && (newForePixel[1] < foreTileHeightBack)){
		
		int x = (int)round(newForePixel[0]);
		int y = (int)round(newForePixel[1]);

		if ((tileDataType == GDT_Byte) && (foreTileDataType == GDT_Byte) &&
		    (foreTileData8UC[y*foreTileWidthBack+x] != foreNoDataVal)){
		  tileData8UC[j*backTileWidthBack+i] = foreTileData8UC[y*foreTileWidthBack+x] - (int)floor(elevationOffset);
		}
		else if ((tileDataType == GDT_Byte) && (foreTileDataType == GDT_Float32) &&
		     !IsFloatNoData(foreTileData32F[y*foreTileWidthBack+x], foreNoDataVal)){
		  tileData8UC[j*backTileWidthBack+i] = (int)floor(255*foreTileData32F[y*foreTileWidthBack+x]);
		}
		else if ((tileDataType == GDT_Float32) &&
			 !IsFloatNoData(foreTileData32F[y*foreTileWidthBack+x], foreNoDataVal)){
		  tileData32F[j*backTileWidthBack+i] = foreTileData32F[y*foreTileWidthBack+x] - elevationOffset;
		}
	      }
	      
	      newForePixel[0]=newForePixel[0]+1;     
	    }
	    newForePixel[1]=newForePixel[1]+1;
	  }
          //fill in the foreground values - END
	  
	  //clear memory
	  GDALClose( (GDALDatasetH)hForeTileDS);

	  //prepare the return params
	  if (tileDataType==GDT_Float32){
            geoImage.f32=tileData32F;
	    geoImage.uc8=NULL;
	    geoImage.dataset = backTileDS;

	  }
	  
	  if (tileDataType==GDT_Byte){
	    geoImage.uc8=tileData8UC;
	    geoImage.f32=NULL;
	    geoImage.dataset = backTileDS;
	  }

	  if (foreTileData32F!=NULL){
	    delete foreTileData32F;
	  }
	  if (foreTileData8UC!=NULL){
	    delete foreTileData8UC;
	  }
	  
	}
      }
    }
    //#endif

    return geoImage;
  }
  //tileParams is given in the original size background image
  //upsampleRatioBackImag is the upsampling factor of the original background
  //image to match the asembled image
  //upsampleRatioForeImg is the upsampling factor of the original foreground
  //image to match the asembled image
  //offset is the elevation offset
  std::vector<std::vector<string> > writeAssembledTiles(GDALDataset* fore, GDALDataset* back, GDALDatasetH hReferenceDS,
							float offsetDEM, int weightingMode,float backUpsamplingFactor,
							int overwriteFlag, int reference,
							std::vector<std::vector<struct TilingParams> > pyrTileParamsArray,
							string mode, string thisResDir)
  {

    cout<<"assembler_util.cc:: writeAssembledTiles(): ASSEMBLER mode= "<<mode<<endl;
    cout<<"assembler_util.cc:: writeAssembledTiles(): backUpsamplingFactor= "<<backUpsamplingFactor<<endl;
    cout<<"assembler_util.cc:: writeAssembledTiles(): overwriteFlag="<<overwriteFlag<<endl;

    int numForeBands = fore->GetRasterCount();
    cout<<"assembler_util.cc:: writeAssembledTiles(): numForeBands="<<numForeBands<<endl;
    
    std::vector<string> assembledTileFilenameArray;
    assembledTileFilenameArray.resize(0);

    vector< vector<string> > filenameArray;
    
    //write tiles only at native resolution of the foreground which is the highest resolution
    //in the pyramid and corresponding to level 0.
    for (unsigned int i= 0; i < pyrTileParamsArray[0].size(); i++){

      struct bbox tilePixelBBoxRef;
       tilePixelBBoxRef.xl = pyrTileParamsArray[0][i].back_xl;
       tilePixelBBoxRef.yt = pyrTileParamsArray[0][i].back_yt;
       tilePixelBBoxRef.xr = pyrTileParamsArray[0][i].back_xr;
       tilePixelBBoxRef.yb = pyrTileParamsArray[0][i].back_yb;
       
       cout<<"assembler_util.cc:: writeAssembledTiles(): tile"
	   <<": "<<pyrTileParamsArray[0][i].horTileIndex
	   <<", "<<pyrTileParamsArray[0][i].verTileIndex
	   <<", xl="<<tilePixelBBoxRef.xl
	   <<", yt="<<tilePixelBBoxRef.yt 
	   <<", xr="<<tilePixelBBoxRef.xr 
	   <<", yb="<<tilePixelBBoxRef.yb <<endl;


       vector<string> filenames;

       //increase the number of bands if there are already multi band background tiles - START
       int numBackBands = 0;
       for (int backBand = 0; backBand < 3; backBand++){
	 stringstream ss;
	 ss<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_"<<backBand+1;
	 string assembledTileFilename;
	 if (mode.compare("DRG")==0){
	   assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_drg.tif";
	 }
	 if (mode.compare("DEM")==0){
	   assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_dem.tif";
	 }
	 //cout<<"assembledTileFilename="<<assembledTileFilename<<endl;
	 if (is_file_exist(assembledTileFilename.c_str())){ 
	   numBackBands++;
	 }
       }
       if (numBackBands==0){
	 numBackBands=1;
       }
       //cout<<"assembler.cc:: writeAssembledTiles(): numBackBands="<<numBackBands<<", "<<mode<<endl;
       //cout<<"assembler.cc:: writeAssembledTiles(): numForeBands="<<numForeBands<<", "<<mode<<endl;
       //increase the number of bands if there are already multi band background tiles - END
       

       if ((numForeBands==1) && (numBackBands==1)){//heavily tested

	 struct GeoImage geoImage;

	 stringstream ss;
	 ss<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_"<<"1";
         string assembledTileFilename;

	 if (mode.compare("DEM")==0){
	     assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_dem.tif";
	 }
	 if (mode.compare("DRG")==0){
	   assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_drg.tif";
	 }
	 
	 geoImage = writeAssembledTile(fore, back, hReferenceDS, offsetDEM,
				       weightingMode, backUpsamplingFactor,
				       tilePixelBBoxRef, mode, 1, reference,
				       assembledTileFilename);
	 
	 if (mode.compare("DEM")==0){
	   if (geoImage.f32!=NULL){
	     saveGDALDataset(geoImage.dataset, geoImage.f32);
	     filenames.push_back(assembledTileFilename);
	     assembledTileFilenameArray.push_back(assembledTileFilename);
	     delete geoImage.f32;
	     //GDALClose( (GDALDatasetH) geoImage.dataset);
	   }
	 }
	 
	 if (mode.compare("DRG")==0){  
	   if (geoImage.uc8!=NULL){
	     saveGDALDataset(geoImage.dataset, geoImage.uc8);
	     filenames.push_back(assembledTileFilename);
	     assembledTileFilenameArray.push_back(assembledTileFilename);
	     delete geoImage.uc8;
	     //GDALClose( (GDALDatasetH) geoImage.dataset);
	   }
	 }
	 
	 GDALClose( (GDALDatasetH) geoImage.dataset);
       }
    

       if ((numForeBands==1) && (numBackBands==3)){//has not been tested yet

	  for (int backBand = 0; backBand<numBackBands; backBand++){

	    struct GeoImage geoImage;
	    
	    stringstream ss;
	    ss<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_"<<backBand+1;
	    string assembledTileFilename;
	    
	    if (mode.compare("DEM")==0){
	      assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_dem.tif";
	    }
	    if (mode.compare("DRG")==0){
	      assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_drg.tif";
	    }
	    
	    geoImage = writeAssembledTile(fore, back, hReferenceDS, offsetDEM,
					  weightingMode, backUpsamplingFactor,
					  tilePixelBBoxRef, mode, 1, reference,
					  assembledTileFilename);

	    if (mode.compare("DEM")==0){
	      if (geoImage.f32!=NULL){
		saveGDALDataset(geoImage.dataset, geoImage.f32);
		filenames.push_back(assembledTileFilename);
		assembledTileFilenameArray.push_back(assembledTileFilename);
		delete geoImage.f32;
		//GDALClose( (GDALDatasetH) geoImage.dataset);
	      }
	    }
	    
	    if (mode.compare("DRG")==0){  
	      if (geoImage.uc8!=NULL){
		saveGDALDataset(geoImage.dataset, geoImage.uc8);
		filenames.push_back(assembledTileFilename);
		assembledTileFilenameArray.push_back(assembledTileFilename);
		delete geoImage.uc8;
		//GDALClose( (GDALDatasetH) geoImage.dataset);
	      }
	    }
	    
	    GDALClose( (GDALDatasetH) geoImage.dataset);
	  }
       }

       if ((numForeBands==3) && (numBackBands==1)){//tested
	 
	 //copy gray back to color bands
	 for (int foreBand = 0; foreBand<numForeBands; foreBand++){

	   struct GeoImage geoImage;
		    
	   string assembledTileFilename;
	   
	   stringstream ss;
	   ss<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_"<<foreBand+1;
	   if (mode.compare("DEM")==0){
	     assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_dem.tif";
	   }
	   if (mode.compare("DRG")==0){
	     assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_drg.tif";
	   }

	   //if color background does not exits, try to make a copy of the gray background.
	   if (!is_file_exist(assembledTileFilename.c_str())){
	     stringstream ss1;
	     ss1<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_1";
	     string assembledTileFilenameGray = thisResDir + string("/assembled_")+ss1.str()+"_drg.tif";
	     string command;
	     command = "cp "+assembledTileFilenameGray+" "+assembledTileFilename;
	     system(command.c_str());
	   }
	   
	   geoImage = writeAssembledTile(fore, back, hReferenceDS, offsetDEM,
					 weightingMode, backUpsamplingFactor,
					 tilePixelBBoxRef, mode, foreBand+1, reference,
					 assembledTileFilename);

	   if (mode.compare("DEM")==0){
	     if (geoImage.f32!=NULL){
	       saveGDALDataset(geoImage.dataset, geoImage.f32);
	       filenames.push_back(assembledTileFilename);
	       assembledTileFilenameArray.push_back(assembledTileFilename);
	       delete geoImage.f32;
	       //GDALClose( (GDALDatasetH) geoImage.dataset);
	     }
	   }
	   
	   if (mode.compare("DRG")==0){  
	     if (geoImage.uc8!=NULL){
	       saveGDALDataset(geoImage.dataset, geoImage.uc8);
	       filenames.push_back(assembledTileFilename);
	       assembledTileFilenameArray.push_back(assembledTileFilename);
	       delete geoImage.uc8;
	       //GDALClose( (GDALDatasetH) geoImage.dataset);
	     }
	     else{
	       //delete the copy of the background file.
	       if (foreBand!=0){
		 string command;
		 command = "rm "+assembledTileFilename;
		 system(command.c_str());
	       }
	     }
	   }
	   
	   GDALClose( (GDALDatasetH) geoImage.dataset);
	 }
       }

       if ((numForeBands==3) && (numBackBands==3)){//tested

	 
	 for (int foreBand = 0; foreBand<numForeBands; foreBand++){

	   struct GeoImage geoImage;
           string assembledTileFilename;		 
	   stringstream ss;
	   ss<<pyrTileParamsArray[0][i].horTileIndex<<"_"<<pyrTileParamsArray[0][i].verTileIndex<<"_"<<foreBand+1;
	   if (mode.compare("DEM")==0){
	     assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_dem.tif";
	   }
	   if (mode.compare("DRG")==0){
	     assembledTileFilename = thisResDir + string("/assembled_")+ss.str()+"_drg.tif";
	   }
	   
	   geoImage = writeAssembledTile(fore, back, hReferenceDS, offsetDEM,
					 weightingMode, backUpsamplingFactor,
					 tilePixelBBoxRef, mode, foreBand+1, reference,
					 assembledTileFilename);
	   
	   if (mode.compare("DEM")==0){
	     if (geoImage.f32!=NULL){
	       saveGDALDataset(geoImage.dataset, geoImage.f32);
	       filenames.push_back(assembledTileFilename);
	       assembledTileFilenameArray.push_back(assembledTileFilename);
	       delete geoImage.f32;
	       //GDALClose( (GDALDatasetH) geoImage.dataset);
	     }
	   }
	   
	   if (mode.compare("DRG")==0){       
	     if (geoImage.uc8!=NULL){
	       saveGDALDataset(geoImage.dataset, geoImage.uc8);
	       filenames.push_back(assembledTileFilename);
	       assembledTileFilenameArray.push_back(assembledTileFilename);
	       delete geoImage.uc8;
	       //GDALClose( (GDALDatasetH) geoImage.dataset);
	     }
	   }
	   
	   GDALClose( (GDALDatasetH) geoImage.dataset);
	 }
       }
      
       
       filenameArray.push_back(filenames);    
    }
    
    return filenameArray;
  }

  
  vector<vector<string> > assembleGeoData(GDALDataset* fore,  GDALDataset* back,
					 struct AssemblerParams assemblerParams, 
					 string mode,
					 struct bbox unionGeoBBox,
					 string foreFilename, string thisResDir,
					 int overwriteFlag)
  {

    ResampleParams resampleParams;
      
    vector<vector<string> > assembledTileFilenameArray;
    std::vector<std::vector<struct TilingParams> > pyrTileParamsArray;
    float backUpsamplingFactor, foreUpsamplingFactor;

    float offsetDEM = 0;

    vector<float> upsampleFactors = computeResampleFactors(fore, back,
							   assemblerParams.backReference, assemblerParams.foreMaxPPD/*,//mpp*/);
    
    //computes background terrain tiling scheme for a given foreground terrain
    if (mode.compare("DEM")==0){
     
      resampleParams.imageROILeft = assemblerParams.roiDEM[0];
      resampleParams.imageROITop  = assemblerParams.roiDEM[1];
      resampleParams.imageROIRight = assemblerParams.roiDEM[2];
      resampleParams.imageROIBottom = assemblerParams.roiDEM[3];
      resampleParams.tileSize = assemblerParams.tileSizeDEM;
      resampleParams.tilePaddingTop = assemblerParams.paddingParamsDEM[0];
      resampleParams.tilePaddingLeft = assemblerParams.paddingParamsDEM[1];
      resampleParams.tilePaddingRight = assemblerParams.paddingParamsDEM[2];
      resampleParams.tilePaddingBottom = assemblerParams.paddingParamsDEM[3];
      resampleParams.pyrResampleFactor=0.5;
      resampleParams.pyramidMode = assemblerParams.pyramidMode;
      
      //determine the offset between the two images;
      offsetDEM = computeOffset(fore, back);
      string regFilename = GetFilenameNoPath(GetFilenameNoExt(foreFilename))+"_reg.txt";
      saveRegistrationParams(offsetDEM, thisResDir+"/"+regFilename);
    }
    
    
    if (mode.compare("DRG")==0){
      
      resampleParams.imageROILeft = assemblerParams.roiDRG[0];
      resampleParams.imageROITop  = assemblerParams.roiDRG[1];
      resampleParams.imageROIRight = assemblerParams.roiDRG[2];
      resampleParams.imageROIBottom = assemblerParams.roiDRG[3];
      resampleParams.tileSize = assemblerParams.tileSizeDRG;
      resampleParams.tilePaddingTop = assemblerParams.paddingParamsDRG[0];
      resampleParams.tilePaddingLeft = assemblerParams.paddingParamsDRG[1];
      resampleParams.tilePaddingRight = assemblerParams.paddingParamsDRG[2];
      resampleParams.tilePaddingBottom = assemblerParams.paddingParamsDRG[3];
      resampleParams.pyrResampleFactor=0.5;
      resampleParams.pyramidMode = assemblerParams.pyramidMode;
    }
 
    pyrTileParamsArray = makeVirtualTiles(back, unionGeoBBox,resampleParams);

    float foreToReferenceResampleFactor;
 
    GDALDatasetH hReferenceDS = makeReferenceLayer(fore, back, unionGeoBBox, assemblerParams,
        					   mode, &foreToReferenceResampleFactor);
    
    assembledTileFilenameArray = writeAssembledTiles(fore, back, hReferenceDS,
						     offsetDEM,  assemblerParams.weightingMode,
						     foreToReferenceResampleFactor, overwriteFlag,
						     assemblerParams.backReference, pyrTileParamsArray,
						     mode, thisResDir);
    
    GDALClose(hReferenceDS);
    
    return assembledTileFilenameArray;
  }

  //assembles DEM or DRG data separately
  vector<vector<string> > assembleGeoData(string foreFilename, string backFilename,
					  struct AssemblerParams assemblerParams, 
					  string mode, string resDir, int overwriteFlag)
  {
    
    string thisResDir;
    std::vector<vector<string> > assembledTileFilenameArray;
    
    //read the fore data  
    GDALDataset* fore = (GDALDataset *) GDALOpen(foreFilename.c_str(), GA_ReadOnly); 
    
    //compute fore bounding box in geographic coordinates
    vector<double> foreBBoxGeographic = computeGeographicBBox(fore);

    vector<string> backTileFilenameArray = AccessDataFilesFromInput(backFilename, "tif");
   
    //select overlaping tiles
    vector<string> overlapBackTileArray;
    overlapBackTileArray = findOverlapFiles(backTileFilenameArray, fore);

    cout<<"numOverlapBackTileArray="<< overlapBackTileArray.size()<<endl;
    vector<float> unionBBoxGeographic;
    unionBBoxGeographic.resize(4);
    unionBBoxGeographic[0]=foreBBoxGeographic[0];
    unionBBoxGeographic[1]=foreBBoxGeographic[1];
    unionBBoxGeographic[2]=foreBBoxGeographic[2];
    unionBBoxGeographic[3]=foreBBoxGeographic[3];
    
    for (unsigned int k = 0; k < overlapBackTileArray.size(); k++){
      
       GDALDataset* tileDS = (GDALDataset *) GDALOpen(overlapBackTileArray[k].c_str(), GA_ReadOnly); 
       vector<double> bbox = computeGeographicBBox(tileDS);
       if (bbox[0]<unionBBoxGeographic[0]){
	 unionBBoxGeographic[0]=bbox[0];
       }
       if (bbox[1]<unionBBoxGeographic[0]){
	 unionBBoxGeographic[1]=bbox[1];
       }
       if (bbox[2]>unionBBoxGeographic[2]){
	 unionBBoxGeographic[2]=bbox[2];
       }
       if (bbox[3]>unionBBoxGeographic[3]){
	 unionBBoxGeographic[3]=bbox[3];
       }
       GDALClose((GDALDatasetH)tileDS);
       cout<<"********************unionBBox "
	   <<unionBBoxGeographic[0]<<", "
	   <<unionBBoxGeographic[1]<<", "
	   <<unionBBoxGeographic[2]<<", "
	   <<unionBBoxGeographic[3]<<endl;
       //}
    struct bbox unionGeoBBox;
    
    /*
    vector<float> unionBBoxGeographic;
    vector<string> DEMFilenameArray;
    DEMFilenameArray.resize(overlapBackDEMFilenameArray.size()+1);
    DEMFilenameArray[0]=foreDEMFilename;
    for (int k=0; k < overlapBackDEMFilenameArray.size(); k++){
      DEMFilenameArray[k+1]=overlapBackDEMFilenameArray[k];
    }
    unionBBoxGeographic = computeUnionGeographicBBox(DEMFilenameArray);
    */
    
    //for (unsigned int k = 0; k < overlapBackTileArray.size(); k++){
      
      //read the background file
      GDALDataset* back = (GDALDataset *) GDALOpen(overlapBackTileArray[k].c_str(), GA_ReadOnly); 
      
      if (mode.compare("DEM")==0){
	//build the result dirname for the assembled tile
	//FIXME: This is a huge hack
	//dirname is rootdirname_modifiedBackfilename
	string tmp=GetFilenameNoPath(GetFilenameNoExt(overlapBackTileArray[k]));            
	FindAndReplace(tmp, "_dem", "");
	FindAndReplace(tmp, "assembled_", "");
	thisResDir = resDir + "_" + tmp;
	
	//create the results directory name
	string command = "mkdir "+ thisResDir;
	cout<<"command="<<command<<endl;
	system(command.c_str());
      }
      
      if (mode.compare("DRG")==0){
	//build the result dirname for the assembled tile
	//FIXME: This is a huge hack
	string tmp=GetFilenameNoPath(GetFilenameNoExt(overlapBackTileArray[k]));            
	FindAndReplace(tmp, "_drg", "");
	FindAndReplace(tmp, "assembled_", "");
	thisResDir = resDir + "_" + tmp;

	//create the results directory name
	string command = "mkdir "+ thisResDir;
	cout<<"command="<<command<<endl;
	system(command.c_str());	
      }
            
      assembledTileFilenameArray = assembleGeoData(fore, back, assemblerParams, mode,
						   unionGeoBBox,
						   foreFilename, thisResDir, overwriteFlag);
      
      GDALClose((GDALDatasetH)back);
    }
    
    GDALClose((GDALDatasetH)fore);
    
    return assembledTileFilenameArray;
  }
  
  vector<struct TerrainFilename >assembleGeoData(GDALDataset *foreDEM, GDALDataset *backDEM,
						 GDALDataset *foreDRG, GDALDataset *backDRG,
						 struct bbox unionGeoBBox,
						 float offsetDEM, struct AssemblerParams assemblerParams, 
						 string thisResDir, int overwriteFlag)
  {
    
    vector<vector<string> >assembledTileFilenameArray;
    assembledTileFilenameArray.resize(2);
    
    ResampleParams demResampleParams;
    demResampleParams.imageROILeft = assemblerParams.roiDEM[0];
    demResampleParams.imageROITop  = assemblerParams.roiDEM[1];
    demResampleParams.imageROIRight = assemblerParams.roiDEM[2];
    demResampleParams.imageROIBottom = assemblerParams.roiDEM[3];
    demResampleParams.tileSize = assemblerParams.tileSizeDEM;
    demResampleParams.tilePaddingTop = assemblerParams.paddingParamsDEM[0];
    demResampleParams.tilePaddingLeft = assemblerParams.paddingParamsDEM[1];
    demResampleParams.tilePaddingRight = assemblerParams.paddingParamsDEM[2];
    demResampleParams.tilePaddingBottom = assemblerParams.paddingParamsDEM[3];
    demResampleParams.pyrResampleFactor=0.5;
    demResampleParams.pyramidMode = assemblerParams.pyramidMode;
    
    ResampleParams drgResampleParams;
    drgResampleParams.imageROILeft = assemblerParams.roiDRG[0];
    drgResampleParams.imageROITop  = assemblerParams.roiDRG[1];
    drgResampleParams.imageROIRight = assemblerParams.roiDRG[2];
    drgResampleParams.imageROIBottom = assemblerParams.roiDRG[3];
    drgResampleParams.tileSize = assemblerParams.tileSizeDRG;
    drgResampleParams.tilePaddingTop = assemblerParams.paddingParamsDRG[0];
    drgResampleParams.tilePaddingLeft = assemblerParams.paddingParamsDRG[1];
    drgResampleParams.tilePaddingRight = assemblerParams.paddingParamsDRG[2];
    drgResampleParams.tilePaddingBottom = assemblerParams.paddingParamsDRG[3];
    drgResampleParams.pyrResampleFactor=0.5;
    drgResampleParams.pyramidMode = assemblerParams.pyramidMode;
  
    vector<float> upsampleFactors;
    float foreToReferenceOffsetDEM;
    float foreToReferenceResampleFactorDEM;

    int numForeBands = foreDRG->GetRasterCount();
    cout<<"assembler_util.cc:: assembleGeoData(): numForeBands="<<numForeBands<<endl;


    numForeBands=foreDRG->GetRasterCount();
    
    //make DEM and DRG data reference layer - END
    
    //*****************************************************************
    //make DEM and DRG virtual reference layer - START
    
    GDALDatasetH hReferenceDEM = makeReferenceLayer(foreDEM, backDEM,
                                                    unionGeoBBox, assemblerParams, "DEM",
                                                    &foreToReferenceResampleFactorDEM);
    cout<<"foreToReferenceResampleFactorDEM="<<foreToReferenceResampleFactorDEM<<endl;
    float foreToReferenceOffsetDRG = 0;
    float foreToReferenceResampleFactorDRG;
    GDALDatasetH hReferenceDRG = makeReferenceLayer(foreDRG, backDRG,
                                                    unionGeoBBox, assemblerParams, "DRG", 
                                                    &foreToReferenceResampleFactorDRG);
    
    cout<<"hReferenceDRG: "<<GDALGetRasterXSize((GDALDataset*)hReferenceDRG)<<", "
	<<GDALGetRasterYSize((GDALDataset*)hReferenceDRG)<<endl;
    cout<<"foreToReferenceResampleFactorDRG="<<foreToReferenceResampleFactorDRG<<endl;
   
    //make DEM and DRG virtual reference layer - END

    //*****************************************************************
    //multi-scale tiling DEM and DRG strategy for the reference layer - START.
    //only bottom layer (0) is used in writeAssembledTiles
    std::vector<std::vector<struct TilingParams> > pyrTileStrategyRefDEM;
    pyrTileStrategyRefDEM = makeVirtualTiles((GDALDataset*) hReferenceDEM,
					     unionGeoBBox,
					     demResampleParams);
    
    std::vector<std::vector<struct TilingParams> > pyrTileStrategyRefDRG;
    pyrTileStrategyRefDRG = changeVirtualTiles((GDALDataset*)hReferenceDEM, demResampleParams,
					       (GDALDataset*)hReferenceDRG, drgResampleParams,
					       pyrTileStrategyRefDEM);
    
    vector<struct TerrainFilename> assembledTerrainFilenameArray;
    
    //*****************************************************************
    //make DEM and DRG data reference layer - START
    if ((assemblerParams.backReference==1)||(assemblerParams.backReference==3)){//back is reference
      //reference DEM uses the sampling grid of background and background resolution
      foreToReferenceOffsetDEM = offsetDEM;
    }
    else{//fore is reference
      foreToReferenceOffsetDEM=0;
      //reproject background on reference (foreground).
      CPLErr err = GDALReprojectImage (backDEM, backDEM->GetProjectionRef(), (GDALDataset*)hReferenceDEM, NULL,
                                       GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL);
      //add DEM offset to the background
      addToGDALDataset((GDALDataset*)hReferenceDEM, -offsetDEM);
    }
    //make DEM and DRG data reference layer - END
 
    //*****************************************************************
    //DEM and DRG tile assembler - START.
    vector<vector<string> > assembledDEMFilenameArray;
    vector<vector<string> > assembledDRGFilenameArray;
    assembledDEMFilenameArray = writeAssembledTiles(foreDEM, backDEM, hReferenceDEM,
						    foreToReferenceOffsetDEM,  assemblerParams.weightingMode, 
						    foreToReferenceResampleFactorDEM, overwriteFlag,
						    assemblerParams.backReference,
						    pyrTileStrategyRefDEM, "DEM", thisResDir);//foreMaxPPD is mpp in fact
    
    assembledDRGFilenameArray = writeAssembledTiles(foreDRG, backDRG, hReferenceDRG,
						    foreToReferenceOffsetDRG,  assemblerParams.weightingMode, 
						    foreToReferenceResampleFactorDRG, overwriteFlag,
						    assemblerParams.backReference,
						    pyrTileStrategyRefDRG, "DRG", thisResDir);

    //DEM and DRG tile assembler - END.

    cout<<"numDEMs="<<assembledDEMFilenameArray.size()<<endl;
    cout<<"numDRGs="<<assembledDRGFilenameArray.size()<<endl;
    
 
    for (int k = 0; k < assembledDEMFilenameArray.size(); k++){
      if ((assembledDEMFilenameArray[k].size()>0)&&
	  (assembledDRGFilenameArray[k].size()>0)){
	struct TerrainFilename assembledTerrainFilename;
	
	for (int b=0; b<assembledDEMFilenameArray[k].size(); b++){
	  assembledTerrainFilename.demFilenames.push_back(assembledDEMFilenameArray[k][b]);
	}
	for (int b=0; b<assembledDRGFilenameArray[k].size(); b++){
	  assembledTerrainFilename.drgFilenames.push_back(assembledDRGFilenameArray[k][b]);
	}
	assembledTerrainFilenameArray.push_back(assembledTerrainFilename);
      }
    }
    
    cout<<"------numAssembledTiles="<<assembledTerrainFilenameArray.size()<<endl;
    for (int k=0;k <assembledTerrainFilenameArray.size();k++){
      for (int b=0; b<assembledTerrainFilenameArray[k].demFilenames.size(); b++){
	cout<<"dem="<<assembledTerrainFilenameArray[k].demFilenames[b]<<endl;
      }
       for (int b=0; b<assembledTerrainFilenameArray[k].drgFilenames.size(); b++){
	cout<<"drg="<<assembledTerrainFilenameArray[k].drgFilenames[b]<<endl;
      }
      cout<<"-----------------------------------"<<endl;
    }

    GDALClose(hReferenceDEM);
    GDALClose(hReferenceDRG);
    
    return assembledTerrainFilenameArray;
  }
  

  vector<struct TerrainFilename> assembleGeoData(string foreDEMFilename, string backDEMFilename,
						 string foreDRGFilename, string backDRGFilename,
						 struct AssemblerParams assemblerParams, 
						 string resDir, int overwriteFlag)
  {

    string thisResDir;

    vector<TerrainFilename> assembledTerrainFilenameArray;
    
    //read the fore data  
    GDALDataset* foreDEM = (GDALDataset *) GDALOpen(foreDEMFilename.c_str(), GA_ReadOnly); 
    GDALDataset* initForeDRG = (GDALDataset *) GDALOpen(foreDRGFilename.c_str(), GA_ReadOnly);
    
    //upsample foreDRG if needed;   
    GDALDatasetH hForeDRG;
    GDALDataset* foreDRG;
    cout<<"resampleFactor="<<assemblerParams.resampleFactorDRG[0]<<endl;

    if (assemblerParams.resampleFactorDRG[0]!=1){
      hForeDRG = resampleGDALDataset(initForeDRG, assemblerParams.resampleFactorDRG[0], "foreDRG.tiff",  string("round"), 0);
      CPLErr err_1 = GDALReprojectImage (initForeDRG, initForeDRG->GetProjectionRef(), hForeDRG, NULL,
					 GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL);
      foreDRG = (GDALDataset*)hForeDRG;
      
      GDALClose((GDALDatasetH)initForeDRG);
      initForeDRG = NULL;
    }
    else{
      foreDRG = initForeDRG;
      hForeDRG = (GDALDatasetH)initForeDRG; 
    }
       
    //compute fore bounding box in geographic coordinates
    vector<double> foreBBoxGeographic = computeGeographicBBox(foreDEM);
   
    //select overlaping tiles
    vector<string> backDEMFilenameArray = AccessDataFilesFromInput(backDEMFilename, "tif"); 
    vector<string> overlapBackDEMFilenameArray;
    overlapBackDEMFilenameArray = findOverlapFiles(backDEMFilenameArray, foreDEM);
    
    vector<string> overlapBackDRGFilenameArray;
    vector<string> backDRGFilenameArray = AccessDataFilesFromInput(backDRGFilename, "tif"); 
    overlapBackDRGFilenameArray = findOverlapFiles(backDRGFilenameArray, foreDRG);

    for (unsigned int k = 0; k < overlapBackDEMFilenameArray.size(); k++){
      
      //read the background file
      GDALDataset* backDEM = (GDALDataset *) GDALOpen(overlapBackDEMFilenameArray[k].c_str(), GA_ReadOnly); 
      GDALDataset* backDRG = (GDALDataset *) GDALOpen(overlapBackDRGFilenameArray[k].c_str(), GA_ReadOnly); 

      
      //compute the union bounding box in geographc coordinates - START
   
      struct bbox unionGeoBBox;
      if (assemblerParams.pyramidMode==1){
         vector<float> unionBBoxGeographic;
	 vector<string> DEMFilenameArray;
	 DEMFilenameArray.resize(2);
	 DEMFilenameArray[0]=foreDEMFilename;
	 DEMFilenameArray[1]=overlapBackDEMFilenameArray[k];
	 
	 unionBBoxGeographic = computeUnionGeographicBBox(DEMFilenameArray);
	 
	 unionGeoBBox.xl = unionBBoxGeographic[0];
	 unionGeoBBox.yb = unionBBoxGeographic[2];
	 unionGeoBBox.xr = unionBBoxGeographic[1];
	 unionGeoBBox.yt = unionBBoxGeographic[3];
      }
      if (assemblerParams.pyramidMode==0){
          vector<float> unionBBoxGeographic;
	  vector<string> DEMFilenameArray;
	  DEMFilenameArray.resize(1);
	  DEMFilenameArray[0]=overlapBackDEMFilenameArray[k];
	  unionBBoxGeographic = computeUnionGeographicBBox(DEMFilenameArray);
	  unionGeoBBox.xl = unionBBoxGeographic[0];
	  unionGeoBBox.yb = unionBBoxGeographic[2];
	  unionGeoBBox.xr = unionBBoxGeographic[1];
	  unionGeoBBox.yt = unionBBoxGeographic[3];
      }
      //compute the union bounding box in geographc coordinates - END
      
      //build the result dirname for the assembled tile
      //FIXME: This is a huge hack
      //dirname is rootdirname_modifiedBackfilename
      string tmp=GetFilenameNoPath(GetFilenameNoExt(overlapBackDEMFilenameArray[k]));            
      FindAndReplace(tmp, "_dem", "");
      FindAndReplace(tmp, "assembled_", "");
      thisResDir = resDir + "_" + tmp;
      cout<<"thisResDir="<<thisResDir<<endl;
      //create the results directory name
      string command = "mkdir "+ thisResDir;
      cout<<"command="<<command<<endl;
      system(command.c_str());

      //registration
      float offsetDEM = 0;
      offsetDEM = computeOffset(foreDEM, backDEM);
      string regFilename = GetFilenameNoPath(GetFilenameNoExt(foreDEMFilename))+"_"
	+ GetFilenameNoPath(GetFilenameNoExt(overlapBackDEMFilenameArray[k]))+"_reg.txt";
      saveRegistrationParams(offsetDEM, thisResDir+"/"+regFilename);
 
      //assembledTerrainFilenameArray
      vector<struct TerrainFilename> assembledTerrainFilenames
	= assembleGeoData(foreDEM, backDEM, foreDRG, backDRG,
						      unionGeoBBox, offsetDEM,
						      assemblerParams, thisResDir, overwriteFlag);

      for (int k = 0; k < assembledTerrainFilenames.size(); k++){
         assembledTerrainFilenameArray.push_back(assembledTerrainFilenames[k]);
      }
      
      //save assembledFilenames to file
      string assembledDEMTileFilenameList=thisResDir+"/"+"dem_list.txt";
      ofstream foutdem(assembledDEMTileFilenameList.c_str());
      for (int i=0; i<assembledTerrainFilenameArray.size(); i++){
          foutdem << assembledTerrainFilenameArray[i].demFilenames[0] << endl;
      }
      foutdem.close();
      

      string assembledDRGTileFilenameList=thisResDir+"/"+"drg_list.txt";
      ofstream foutdrg(assembledDRGTileFilenameList.c_str());
      for (int i=0; i<assembledTerrainFilenameArray.size(); i++){
	for (int b=0; b<assembledTerrainFilenameArray[i].drgFilenames.size(); b++){
          foutdrg << assembledTerrainFilenameArray[i].drgFilenames[b] << endl;
	}
      }
      foutdrg.close();
 
      GDALClose((GDALDatasetH)backDEM);
      GDALClose((GDALDatasetH)backDRG);

    }

    GDALClose((GDALDatasetH)foreDEM);

    if (initForeDRG!=NULL){
        GDALClose((GDALDatasetH)hForeDRG);
    }
    return assembledTerrainFilenameArray;
  }




  
#if 0
  void addBand(string resDirname, struct TilingParams &tileParams, string fileType, int bandIndex)
  {
    stringstream ssInput;
    ssInput<<tileParams.horTileIndex<<"_"<<tileParams.verTileIndex<<"_"<<fileType;
    string outputFilename = resDirname + string("/assembled_")+ssInput.str()+".tif";
    
    stringstream ssOutput;
    ssOutput<<tileParams.horTileIndex<<"_"<<tileParams.verTileIndex<<"_"<<bandIndex<<"_"<<fileType;
    string inputFilename = resDirname + string("/assembled_")+ssOutput.str()+".tif";

    //read band tile
    cout<<"inputFilename="<<inputFilename<<endl;
    GDALDataset* inputDS = (GDALDataset *) GDALOpen(inputFilename.c_str(), GA_ReadOnly);
    int inputWidth = GDALGetRasterXSize(inputDS);
    int inputHeight = GDALGetRasterYSize(inputDS);
    GDALRasterBand* inputBand = inputDS->GetRasterBand(1);
    GDALDataType inputType = inputBand->GetRasterDataType();

    float *inputData32F = NULL;
    unsigned char *inputData8UC = NULL;

    if (inputType == GDT_Byte){ 
      inputData8UC = new unsigned char[inputHeight*inputWidth];
      inputBand->RasterIO(GF_Read, 0, 0, 
			  inputWidth, inputHeight,  
			  inputData8UC, inputWidth, inputHeight, 
			  inputType, 0, 0);
      //foreNoDataVal = 0.0;
      //foreDSNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)foreBand, &pbSuccess);
    }
    else if (inputType == GDT_Float32) {
      inputData32F = new float[inputHeight*inputWidth];
      inputBand->RasterIO(GF_Read, 0, 0, 
			  inputWidth, inputHeight,  
			  inputData32F, inputWidth, inputHeight, 
			  inputType, 0, 0);
      //foreNoDataVal = NAN;
      //foreDSNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)foreBand, &pbSuccess);
    }

    GDALDataset* outputDS = (GDALDataset *) GDALOpen(outputFilename.c_str(), GA_Update);

    if (outputDS==NULL){
      GDALDriverH hOutputDriver = GDALGetDriverByName("GTiff");
      CPLErr imError = GDALRenameDataset (hOutputDriver, outputFilename.c_str(),
					  inputFilename.c_str());	
    } 
    else{
      cout<<"outputFilename="<<outputFilename<<endl;

      //create a virtual dataset from outputDS
      GDALDriverH hVrtDriver = GDALGetDriverByName("vrt");
      GDALDatasetH outputVrtDSH = GDALCreateCopy(hVrtDriver, /*outputFilename.c_str()*/"out.vrt", (GDALDatasetH)outputDS, FALSE, NULL, NULL, NULL);    
      GDALDataset* outputVrtDS = (GDALDataset*)outputVrtDSH; 
      
      //add a band
      outputVrtDS->AddBand(inputType); 
      GDALRasterBandH hVrtDSBandout = GDALGetRasterBand(outputVrtDSH, bandIndex);

      //TODO: copy the virtual dataset to a geotif dataset
      GDALDriverH  hOutputDriver = GDALGetDriverByName("GTiff");
      GDALDatasetH outputDSH = GDALCreateCopy(hOutputDriver, outputFilename.c_str(), (GDALDatasetH)outputVrtDS, FALSE, NULL, NULL, NULL); 
      GDALDataset* outputDS = (GDALDataset*)outputDSH; 

      GDALClose( (GDALDatasetH) outputVrtDS);

      //populate and save the new band
      GDALRasterBandH hBandout = GDALGetRasterBand(outputDS, bandIndex);
      if (bandIndex==1){
	GDALSetRasterColorInterpretation(hBandout, GCI_BlueBand);
      } 
      if (bandIndex==2){
	GDALSetRasterColorInterpretation(hBandout, GCI_GreenBand);
      } 
      if (bandIndex==3){
	GDALSetRasterColorInterpretation(hBandout, GCI_RedBand);
      } 			
      CPLErr imError = CE_None;
      if (inputType==GDT_Float32){
	imError = GDALRasterIO( hBandout, GF_Write, 0, 0, inputWidth, inputHeight,
				inputData32F, inputWidth , inputHeight , inputType, 0,0);
      }
      if (inputType==GDT_Byte){
	imError = GDALRasterIO( hBandout, GF_Write, 0, 0, inputWidth, inputHeight,
				inputData8UC, inputWidth , inputHeight , inputType, 0,0);
      }
      
      GDALClose( (GDALDatasetH) outputDS);
      
    }
    //save info to the combined file - END
   
    //dealocate memory
    if (inputData32F!=NULL){
      delete inputData32F;
    }
    if (inputData8UC!=NULL){
      delete inputData8UC;
    }

    GDALClose( (GDALDatasetH) inputDS);

  }
#endif
  
#if 0
    float ComputeDEMAccuracy(GeoReference foreGeo, Vector2 forePix, float backAccuracy)
    {
   
        Vector2 forePoint;
        forePoint = foreGeo.pixel_to_point(forePix);
  
        float distToCam = sqrt(forePoint(0)*forePoint(0) + forePoint(1)*forePoint(1));
        float foreWeight = 0.5;
  
        //delta_d = 43 micro, f = 43mm, b = 0.3m
        float foreAccuracy = distToCam*distToCam/300.0; ///delta_r = r*r*delta_d/(b*f)
        foreWeight = backAccuracy/(foreAccuracy+backAccuracy);
  
        return foreWeight;

    }
#endif




    //============================================================================
#if 0
    void SaveAssembledPC(string DEMFilename, string assembledPCFilename)
    {
        cout<<"Saving the PC to file ..."<<endl;
  
        DiskImageView<float>   DEM(DEMFilename);
        GeoReference Geo;
        read_georeference(Geo, DEMFilename);
        int verStep = 20;
        int horStep = 20;
        float noDataVal = -37687.0;
  
 
        ImageViewRef<float>interpDEM = interpolate(edge_extend(DEM.impl(),
                                                               ConstantEdgeExtension()),
                                                   BilinearInterpolation());

        Vector2 pix_0(0,0);
        Vector2 lonlat_0= Geo.pixel_to_lonlat(pix_0);
        Vector3 lonlat3_0(lonlat_0(0), lonlat_0(1), (interpDEM.impl())(0,0));
        Vector3 xyz_0 = Geo.datum().geodetic_to_cartesian(lonlat3_0);
  
        FILE* fp = fopen(assembledPCFilename.c_str(), "wt");
   
        for (int j = 0; j < DEM.impl().rows(); j=j+verStep)
	{
            for (int i = 0; i < DEM.impl().cols(); i=i+horStep)
	    {
                if (!isnan(DEM.impl()(i, j)) &&
		    (DEM.impl()(i, j)) != -noDataVal)
		{
                    Vector2 pix(i,j);
                    Vector2 lonlat = Geo.pixel_to_lonlat(pix);
                    Vector3 lonlat3(lonlat(0), lonlat(1), (interpDEM.impl())(i,j));
                    Vector3 xyz = Geo.datum().geodetic_to_cartesian(lonlat3);
                    fprintf(fp, "%f %f %f\n", xyz[0]-xyz_0[0], xyz[1]-xyz_0[1], xyz[2]-xyz_0[2]);
                }
            }
        }

        fclose(fp);
    }

    //determines the new point after translation by xyz (z is not used)
    //assumes transaltion is x towards North, y towards East, z towards the center of the planet
    Vector2 translateDEMOriginPoint(Vector2 forePointOrig, Vector3 foreXYZOffset, GeoReference const &Geo)
    {
        Vector2 newOrigPoint;
    
        //account for XYZ offset from landing site - START
        cout<<"foreXYZOffset="<<foreXYZOffset<<endl;
        Matrix<double> H_dem;
        H_dem = Geo.transform();
        Vector2 deltaPix;//offset in pixels
        deltaPix(0) = foreXYZOffset(1)/H_dem(0,0);
        deltaPix(1) = foreXYZOffset(0)/H_dem(1,1);
        Vector2 origPix = Geo.point_to_pixel(forePointOrig);
        Vector2 newPix=origPix+deltaPix;
        newOrigPoint = Geo.pixel_to_point(newPix);    
        cout<<"newPoint="<<newOrigPoint<<endl;
   
        return newOrigPoint;
    }
#endif
}





