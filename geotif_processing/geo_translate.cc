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


#include "gdalwarper.h"
#include "GeoUtils.h"

using namespace std;
using namespace at;

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) 
{

  std:: string filename = argv[1];
  double tiePointX = atof(argv[2]);
  double tiePointY = atof(argv[3]);
  
  GDALAllRegister();
     
  //read the fore data  
  GDALDataset* src = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
  int width  = GDALGetRasterXSize(src);
  int height = GDALGetRasterYSize(src);
  
  //translate the header
  std::vector<double> newTiePoint;
  newTiePoint.resize(2);
  newTiePoint[0] = tiePointX;
  newTiePoint[1] = tiePointY;
  
  GDALDatasetH dstH = translateGDALDataset(src, newTiePoint, "translated.tiff");
  CPLErr err = GDALReprojectImage (src, src->GetProjectionRef(), dstH, NULL,
  				   GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL);
  CPLErr imError   = CE_None;

  GDALRasterBand* band;
  int numBands = src->GetRasterCount();
  GDALRasterBandH hBandout;
  
  if (numBands==1){
    band = src->GetRasterBand(1);
    GDALDataType dataType=band->GetRasterDataType();
    
    float* dataVals = new float[width*height];
    band->RasterIO(GF_Read, 0, 0, width, height,
		   dataVals, width, height, 
		   dataType, 0, 0);

    hBandout= GDALGetRasterBand( dstH, 1 );
    imError = GDALRasterIO( hBandout, GF_Write, 0, 0, width, height,
			    dataVals, width , height , dataType ,
			    0, 0);
    delete dataVals;
  }
  
  if (numBands==3){

    for (int i = 0; i < 3; i++){
      
      band = src->GetRasterBand(i+1);
      GDALDataType dataType=band->GetRasterDataType();
      unsigned char* dataVals = new unsigned char[width*height];
      
      band->RasterIO(GF_Read, 0, 0, width, height,
		     dataVals, width, height, 
		     dataType, 0, 0);
      
      hBandout= GDALGetRasterBand(dstH, i+1);
      imError = GDALRasterIO( hBandout, GF_Write, 0, 0, width, height,
			      dataVals, width, height , dataType ,
			      0, 0);

      delete dataVals;
      
    }

  }
  
  GDALClose(dstH);
  GDALClose((GDALDatasetH)src);
    
  return 0;
}

















