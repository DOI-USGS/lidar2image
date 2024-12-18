// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "gdal_priv.h"
#include "gdalwarper.h"

#include "../common/StringUtils.h"
#include "../geotif_processing/GeoUtils.h"

using namespace std;
using namespace at;

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) 
{
  string filename;
  string resDirname = "";
  float rescaleFactor;
  int filenameFound = 0;
  int rescaleFactorFound = 0;
  int resDirnameFound = 0;
  int testFlag = 0;
  int testError = 0;

  for (int i = 1; i < argc; i++) { //We will iterate over argv[] to get the parameters stored inside.
    //Note that we're starting on 1 because we don't need to know the
    //path of the program, which is stored in argv[0]
    if (i + 1 != argc) { // Check that we haven't finished parsing already
      
      if (string(argv[i]).compare("-t")==0) {
	testFlag = 1;
      }
      else if (string(argv[i]).compare("-i")==0) {
	filename = argv[i + 1];
	filenameFound = 1;
      }
      else if (string(argv[i]).compare("-f")==0) {
	rescaleFactor = atof(argv[i + 1]);
	rescaleFactorFound = 1;
      }
      else if (string(argv[i]).compare("-r")==0) {
	resDirname = argv[i + 1];
	resDirnameFound = 1; 
      }
    }
  }
  
  int usageError = 0;
 
  if (filenameFound==0){
    cout<<"Error! Input filename not found"<<endl;
    usageError = 1;
  }
  if (rescaleFactorFound==0){
    cout<<"Error! rescale factor not found"<<endl;
    usageError = 1;
  }
   
  if (usageError == 1){
    std::cout << "Usage is "<< argv[0]<< " -i inputFilename -f rescaleFactor"
	      <<" [-t] [-r resultDirname]"<<endl;
    std::cin.get();
    exit(0);
  }
  string outFilename = resDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(filename))+"_rescaled.tif";
  
  GDALAllRegister();
     
  //read the data  
  GDALDataset* data = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly); 

  int width  = GDALGetRasterXSize(data);
  int height = GDALGetRasterYSize(data);
  GDALRasterBand* bandData = data->GetRasterBand(1);
  GDALDataType typeData    = bandData->GetRasterDataType();

  float * f32Data = NULL; 
  unsigned char *uC8Data = NULL;
  double noDataVal;
  
  if (typeData == GDT_Byte){ 
    uC8Data = new unsigned char[height*width];
    bandData->RasterIO(GF_Read, 0, 0, 
		       width, height,  
		       uC8Data, width, height, 
		       typeData, 0, 0);
    int pbSuccess = 0;
    noDataVal = 0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData, &pbSuccess);
    if (pbSuccess == 1){
      noDataVal= dSNoDataVal;
    } 	
  }

  if (typeData == GDT_Float32){ 
    f32Data = new float[height*width];
    bandData->RasterIO(GF_Read, 0, 0, 
		       width, height,  
		       f32Data, width, height, 
		       typeData, 0, 0);
    int pbSuccess = 0;
    noDataVal = 0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData, &pbSuccess);
    if (pbSuccess == 1){
      noDataVal = dSNoDataVal;
    } 	
  }
  
  GDALDatasetH resampledDataH = resampleGDALDataset(data, rescaleFactor, outFilename, string("round"), 1);
  GDALDataset* resampledData=(GDALDataset*)resampledDataH;
  GDALRasterBandH hBandout = GDALGetRasterBand(resampledData, 1);
  int resampledDataWidth  = GDALGetRasterXSize(resampledData);
  int resampledDataHeight = GDALGetRasterYSize(resampledData);
  CPLErr imError   = CE_None;
   
  if (typeData == GDT_Byte){ 
    unsigned char*resampled = new unsigned char[resampledDataWidth * resampledDataHeight];  
    bandData->RasterIO(GF_Read, 0, 0, width, height, resampled,
		       resampledDataWidth, resampledDataHeight, GDT_Byte, 0, 0);
    imError = GDALRasterIO( hBandout, GF_Write, 0, 0, resampledDataWidth, resampledDataHeight,
			    resampled, resampledDataWidth , resampledDataHeight , GDT_Byte, 0, 0);
  }
  
  if (typeData == GDT_Float32){
    float*resampled = new float[resampledDataWidth * resampledDataHeight];
    bandData->RasterIO(GF_Read, 0, 0, width, height, resampled,
		       resampledDataWidth, resampledDataHeight, GDT_Float32, 0, 0);
    imError = GDALRasterIO(hBandout, GF_Write, 0, 0, resampledDataWidth, resampledDataHeight,
			   resampled, resampledDataWidth , resampledDataHeight , GDT_Float32, 0, 0);
  }
  
  GDALClose(resampledData);
  GDALClose(data);

 
  return 0;
}

















