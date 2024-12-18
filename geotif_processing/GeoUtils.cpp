// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//#include <errno.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "gdal_priv.h"
#include "gdalwarper.h"
#include <ogr_spatialref.h>

#include "../common/Tiling.h"
#include "GeoUtils.h"

using namespace std;

namespace at 
{

  //upsample/downsample GDALDataset
  GDALDatasetH resampleGDALDataset(GDALDataset* src, float resampleFactor,
				   string dstFilename, string resizeType, int saveFileFlag)
  {

    GDALDatasetH dstH;

    GDALDriverH hDriver;  
    char **pzOptions = NULL;
    int imageBpp = 1;
    
    //create the upsampledData
    hDriver = GDALGetDriverByName("GTiff");
    int resampledWidth, resampledHeight;
    if (resizeType.compare("round")==0){
      resampledWidth=(int)round(GDALGetRasterXSize(src)*resampleFactor);
      resampledHeight=(int)round(GDALGetRasterYSize(src)*resampleFactor);
    }
    if (resizeType.compare("floor")==0){
      resampledWidth=(int)floor(GDALGetRasterXSize(src)*resampleFactor);
      resampledHeight=(int)floor(GDALGetRasterYSize(src)*resampleFactor);
    }
    if (resizeType.compare("ceil")==0){
      resampledWidth=(int)ceil(GDALGetRasterXSize(src)*resampleFactor);
      resampledHeight=(int)ceil(GDALGetRasterYSize(src)*resampleFactor);
    }
    cout<<"resampled size="<<resampledWidth<<", "<<resampledHeight<<endl;
    GDALRasterBand* srcBand=src->GetRasterBand(1);
    GDALDataType resampledDataType=srcBand->GetRasterDataType();
    
    dstH = GDALCreate(hDriver, dstFilename.c_str(),
		      resampledWidth, resampledHeight,
		      imageBpp, resampledDataType, pzOptions);

    cout<<"GeoUtils.cc:: resampledGDALDataset(): resampled size="<<resampledWidth<<", "<<resampledHeight<<endl;
    //copy projection from backData
    GDALSetProjection(dstH, src->GetProjectionRef());
    
    //set geoTransform
    double srcGeoTransform[6];
    src->GetGeoTransform(srcGeoTransform);
        
    double resampledGeoTransform[6];
    vector<double> projCoordTL;
    vector<double> pixelTL; pixelTL.resize(2);
    pixelTL[0]=0;
    pixelTL[1]=0;
    projCoordTL=pixelToProjected(pixelTL, src);
    pixelTL.resize(0);
    
    resampledGeoTransform[0]=projCoordTL[0];
    resampledGeoTransform[1]=srcGeoTransform[1]/resampleFactor;
    resampledGeoTransform[2]=srcGeoTransform[2];
    resampledGeoTransform[3]=projCoordTL[1];
    resampledGeoTransform[4]=srcGeoTransform[4];
    resampledGeoTransform[5]=srcGeoTransform[5]/resampleFactor;

    GDALSetGeoTransform(dstH, resampledGeoTransform );

    if (saveFileFlag == 0){
      string command="rm "+dstFilename;
      system(command.c_str());
    }
    
    return dstH;
  }

  float* uint16ToF32(GDALDataset *data)
  {
    GDALDataset *data32F;
    GDALDatasetH hData32F;
    GDALRasterBandH hBandData32F;
    
    int width  = GDALGetRasterXSize(data);
    int height = GDALGetRasterYSize(data);

    cout<<"width="<<width<<", height="<<height<<endl;
    
    GDALRasterBand* bandData = data->GetRasterBand(1);
    GDALDataType typeData  = bandData->GetRasterDataType();
    short int *imageC16 = new short int[height*width];
    bandData->RasterIO(GF_Read, 0, 0, width, height,  
			  imageC16, width, height, 
			  typeData, 0, 0);

    float noDataValue = bandData->GetNoDataValue();
   
    signed short maxVal = bandData->GetMaximum();
    signed short minVal = bandData->GetMinimum(); 
  
    cout<<"minVal="<<minVal<<", maxVal="<<maxVal<<endl;
    
    double gain = 1.0/(maxVal-minVal);
    double bias = -gain*minVal;
    cout<<"GeoUtils.cc:: Short Int gain="<<gain<<", bias="<<bias<<endl;

    //create a float gdal image - START
    float* imageF32 = new float[width*height];
    for (int i = 0; i < height; i++){
      for (int j = 0; j < width; j++){
	if (imageC16[i*width+j]==noDataValue){
	    imageF32[i*width+j] =  -3.40282265508890445e+38;  
	}
	else{
	    imageF32[i*width+j] = gain*imageC16[i*width+j]+bias;
	}
      }
    }
    delete imageC16;

    return imageF32;
  }
  
  GDALDataset* uint16ToF32GDAL(GDALDataset * data, string filename, int saveFileFlag)
  {
    GDALDataset *data32F;
    GDALDatasetH hData32F;
    GDALRasterBandH hBandData32F;
    int width  = GDALGetRasterXSize(data);
    int height = GDALGetRasterYSize(data);
  
    GDALRasterBand* bandData = data->GetRasterBand(1);
    GDALDataType typeData  = bandData->GetRasterDataType();
    short int *imageC16 = new short int[height*width];
    bandData->RasterIO(GF_Read, 0, 0, width, height,  
			  imageC16, width, height, 
			  typeData, 0, 0);

    float noDataValue = bandData->GetNoDataValue();
   
    signed short maxVal = bandData->GetMaximum();
    signed short minVal = bandData->GetMinimum(); 
  
    cout<<"minVal="<<minVal<<", maxVal="<<maxVal<<endl;
    
    double gain = 1.0/(maxVal-minVal);
    double bias = -gain*minVal;
    cout<<"GeoUtils.cc:: Short Int gain="<<gain<<", bias="<<bias<<endl;

    //create a float gdal image - START
    float* imageF32 = new float[width*height];
    for (int i = 0; i < height; i++){
      for (int j = 0; j < width; j++){
	if (imageC16[i*width+j]==noDataValue){
	    imageF32[i*width+j] =  -3.40282265508890445e+38;  
	}
	else{
	  imageF32[i*width+j] = gain*imageC16[i*width+j]+bias;
	}
      }
    }
    delete imageC16;
    
    GDALDriverH hDriver;  
    hDriver = GDALGetDriverByName("GTiff");
    
    char **pzOptions=NULL;
    
    int numBands = data->GetRasterCount();
    
    hData32F = GDALCreate(hDriver, filename.c_str(),				     
			     width, height,
			     numBands, GDT_Float32, pzOptions);
    
    GDALSetProjection( hData32F, data->GetProjectionRef());
    
    double srcGeoTransform[6];
    data->GetGeoTransform(srcGeoTransform);
    GDALSetGeoTransform( hData32F, srcGeoTransform );
    
    data32F = (GDALDataset*) hData32F;
    
    hBandData32F = GDALGetRasterBand(data32F, 1);
    GDALRasterBand* bandData32F =((GDALDataset*)hData32F)->GetRasterBand(1);
    bandData32F->SetNoDataValue(-3.40282265508890445e+38);
    
    bandData32F->RasterIO(GF_Write, 0, 0, width, height,
			     imageF32, width, height, GDT_Float32, 0, 0);


    if (saveFileFlag == 0){
      string command="rm "+filename;
      system(command.c_str());
    }
    
    return data32F;
 }
  //change layer count
  GDALDatasetH changeBandsGDALDataset(GDALDataset* src, int numRasters, string dstFilename)
  {
    GDALDatasetH hDst;
    GDALDriverH hDriver;  
    char **pzOptions = NULL;

    hDriver = GDALGetDriverByName("GTiff");
    int dstWidth=GDALGetRasterXSize(src);
    int dstHeight=GDALGetRasterYSize(src);
    
    GDALRasterBand* srcBand = src->GetRasterBand(1);
    GDALDataType dstDataType = srcBand->GetRasterDataType();
    
    hDst = GDALCreate(hDriver, dstFilename.c_str(),
		      dstWidth, dstHeight,
		      numRasters, dstDataType, pzOptions);

    //copy projection from backData
    GDALSetProjection( hDst, src->GetProjectionRef());

    //copy the geoTransform;
    double srcGeoTransform[6];
    src->GetGeoTransform(srcGeoTransform);
    GDALSetGeoTransform( hDst, srcGeoTransform);
    return hDst;
  }
  
  
  
  //resizes an existing GDALDataset to a new pixel bounding box.
  //replaces cropGDALDataset
  GDALDatasetH resizeGDALDataset(GDALDataset* src, struct bbox pixelBBox, string dstFilename)
  {
   
    GDALDatasetH hDst;
    GDALDriverH hDriver;  
    char **pzOptions = NULL;
    pzOptions = CSLSetNameValue( pzOptions, "BIGTIFF", "YES" );
  
    hDriver = GDALGetDriverByName("GTiff");
    int dstWidth = fabs(pixelBBox.xr-pixelBBox.xl);
    int dstHeight = fabs(pixelBBox.yb-pixelBBox.yt);
    GDALRasterBand* srcBand = src->GetRasterBand(1);
    GDALDataType dstDataType = srcBand->GetRasterDataType();
      
    hDst = GDALCreate(hDriver, dstFilename.c_str(),
		      dstWidth, dstHeight,
		      src->GetRasterCount(), dstDataType, pzOptions);
    
    //copy projection from backData
    GDALSetProjection( hDst, src->GetProjectionRef());
    
    //set geoTransform
    double srcGeoTransform[6];
    src->GetGeoTransform(srcGeoTransform);
    
    vector<double> pixelTL; pixelTL.resize(2);
    pixelTL[0] = pixelBBox.xl;
    pixelTL[1] = pixelBBox.yt;
    vector<double> projCoordTL;
    projCoordTL = pixelToProjected(pixelTL, src);
    
    double dstGeoTransform[6];
    dstGeoTransform[0]=projCoordTL[0];
    dstGeoTransform[1]=srcGeoTransform[1];
    dstGeoTransform[2]=srcGeoTransform[2];
    dstGeoTransform[3]=projCoordTL[1];
    dstGeoTransform[4]=srcGeoTransform[4];
    dstGeoTransform[5]=srcGeoTransform[5];
        
    GDALSetGeoTransform( hDst, dstGeoTransform );

    string command="rm "+dstFilename;
    //system(command.c_str());

    return hDst;
  }


 
  //resize and resample only the header GDALDatasetH
  GDALDatasetH cropGDALDataset(GDALDataset* origData,  struct bbox &pixelBBox,
			       int bandIndex, float upsampleRatio, string filename)
  {
	
    GDALDatasetH hTileDS;

    int tileWidth, tileHeight;
    tileWidth =round((pixelBBox.xr-pixelBBox.xl)*upsampleRatio);
    tileHeight=round((pixelBBox.yb-pixelBBox.yt)*upsampleRatio);
    
    double tileGeoTransform[6];
    GDALDriverH hDriver;  
    char **pzOptions=NULL;
    pzOptions = CSLSetNameValue( pzOptions, "BIGTIFF", "YES" ); 
    GDALRasterBand* band;
    
    int numBands = origData->GetRasterCount();
    band = origData->GetRasterBand(1);

    GDALDataType tileDataType=band->GetRasterDataType();
    
    hDriver = GDALGetDriverByName("GTiff");
    hTileDS  = GDALCreate(hDriver, filename.c_str(),
			  tileWidth, tileHeight,
			  numBands, tileDataType, pzOptions);
    
    //copy projection from backData
    GDALSetProjection( hTileDS, origData->GetProjectionRef());
    
    //set geoTransform - START
    double origGeoTransform[6];
    origData->GetGeoTransform(origGeoTransform);
    
    vector<double> pixelTL;
    pixelTL.resize(2);
    
    pixelTL[0]=pixelBBox.xl;
    pixelTL[1]=pixelBBox.yt;
    
    vector<double> projCoordTL;
    projCoordTL=pixelToProjected(pixelTL, origData);
    
    tileGeoTransform[0]=projCoordTL[0];
    tileGeoTransform[1]=origGeoTransform[1]/upsampleRatio;
    tileGeoTransform[2]=origGeoTransform[2];
    tileGeoTransform[3]=projCoordTL[1];
    tileGeoTransform[4]=origGeoTransform[4];
    tileGeoTransform[5]=origGeoTransform[5]/upsampleRatio;
    
    pixelTL.resize(0);
    
    GDALSetGeoTransform( hTileDS, tileGeoTransform );
    //set geoTransform - END
     
    return hTileDS;
  
  }

  GDALDatasetH translateGDALDataset(GDALDataset* src, vector<double> newTiePoint,
				    string filename)
  {
    GDALDatasetH hDst;


    int width  = GDALGetRasterXSize(src);
    int height = GDALGetRasterYSize(src);
    cout<<"width="<<width<<", height="<<height<<endl;
    char **pzOptions=NULL;
    pzOptions = CSLSetNameValue( pzOptions, "BIGTIFF", "YES" ); 
 
    int numBands = src->GetRasterCount();

    GDALRasterBand* band;
    band = src->GetRasterBand(1);
    GDALDataType dataType=band->GetRasterDataType();
    
    GDALDriverH hDriver;  
    hDriver = GDALGetDriverByName("GTiff");

    hDst = GDALCreate(hDriver, filename.c_str(),
		      width, height,
		      numBands, dataType, pzOptions);
    
    //copy projection from backData
    GDALSetProjection( hDst, src->GetProjectionRef());
    
 
    double srcGeoTransform[6];
    src->GetGeoTransform(srcGeoTransform);

    double dstGeoTransform[6];
    dstGeoTransform[0] = newTiePoint[0];
    dstGeoTransform[1] = srcGeoTransform[1];
    dstGeoTransform[2] = srcGeoTransform[2];
    dstGeoTransform[3] = newTiePoint[1];
    dstGeoTransform[4] = srcGeoTransform[4];
    dstGeoTransform[5] = srcGeoTransform[5];
       
    GDALSetGeoTransform( hDst, dstGeoTransform );
    //set geoTransform - END
     
    return hDst;
  }
  
  void addToGDALDataset(GDALDataset* src, float val)
  {

    int hRefDEMWidth = GDALGetRasterXSize(src);
    int hRefDEMHeight = GDALGetRasterYSize(src);
    
    GDALRasterBandH hRefDEMBand = GDALGetRasterBand(src, 1);
    float *hRefDEMData32F=new float[hRefDEMWidth*hRefDEMHeight];
    CPLErr imError   = CE_None;
    imError = GDALRasterIO( hRefDEMBand, GF_Read, 0, 0, hRefDEMWidth, hRefDEMHeight,
			    hRefDEMData32F, hRefDEMWidth, hRefDEMHeight, GDT_Float32, 0, 0);
    for (int row=0; row<hRefDEMHeight; row++){
       for (int col=0; col<hRefDEMWidth; col++){
	 hRefDEMData32F[row*hRefDEMWidth+col]=hRefDEMData32F[row*hRefDEMWidth+col]+val;
       }
    }

    imError = GDALRasterIO( hRefDEMBand, GF_Write, 0, 0, hRefDEMWidth, hRefDEMHeight,
			    hRefDEMData32F, hRefDEMWidth, hRefDEMHeight, GDT_Float32, 0, 0);
    delete hRefDEMData32F;
  }
  

struct bbox Get2DBoundingBox(GDALDataset *data)
{
  struct bbox bBoxProjected;
  
  double geoTransform[6] = { 999.0, 999.0, 999.0, 999.0, 999.0, 999.0 };
  if (data->GetGeoTransform(geoTransform) == CE_Failure){
    cerr << "WARNING in geotof_utilcc:: computeGeographicBBox(): could not read geoTransform "
	 << "using:\n\t["
	 << geoTransform[0] << ", " << geoTransform[1] << ", "
	 << geoTransform[2] << ", " << geoTransform[3] << ", "
	 << geoTransform[4] << ", " << geoTransform[5] << "]" << endl;
  }

  int width = GDALGetRasterXSize(data);
  int height = GDALGetRasterYSize(data);
  
  bBoxProjected.xl = geoTransform[0];
  bBoxProjected.xr = geoTransform[0] + geoTransform[1]*width;
  bBoxProjected.yt = geoTransform[3];
  bBoxProjected.yb = geoTransform[3] + geoTransform[5]*height;  
  
  return bBoxProjected;
}

  float* GetNoDataValue(GDALDataset* data)
  {
    float noDataValues[3];
    noDataValues[0] = 0;
    noDataValues[1] = 0;
    noDataValues[2] = 0;

    double matchDEMNoDataVal;
    int pbSuccess = 0;
    int numMatchDEMBands = data->GetRasterCount();
    cout<<"numMatchDEMBands="<<numMatchDEMBands<<endl;
    GDALRasterBand* matchDEMBand;
    if (numMatchDEMBands == 1){
      matchDEMBand = data->GetRasterBand(1);
      matchDEMNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)matchDEMBand, &pbSuccess);
      if (pbSuccess){
	noDataValues[0] = matchDEMNoDataVal;
	noDataValues[1] = matchDEMNoDataVal;
	noDataValues[2] = matchDEMNoDataVal;
      }
    }
    else{
      for (int bandIndex = 0; bandIndex<numMatchDEMBands; bandIndex++){
	matchDEMBand = data->GetRasterBand(bandIndex+1);
	matchDEMNoDataVal = GDALGetRasterNoDataValue((GDALRasterBandH)matchDEMBand, &pbSuccess);
	cout<<"matchDEMNoDataVal="<<matchDEMNoDataVal<<", pbSucces="<<pbSuccess<<endl;
	if ((pbSuccess == 1)&&(bandIndex==0)){
	  noDataValues[0] = matchDEMNoDataVal;
	}
	if ((pbSuccess == 1)&&(bandIndex==1)){
	  noDataValues[1] = matchDEMNoDataVal;
	}
	if ((pbSuccess == 1)&&(bandIndex==2)){
	  noDataValues[2] = matchDEMNoDataVal;
	}
      }
    }
    return noDataValues;
  }
  
  //computes the offset between two gdaldatasets. for dems this is
  //elevation offset, and for ortho-phtos this is the intensity offset
  float computeOffset(GDALDataset* foreData, GDALDataset* backData)
  {
    //read the foreground - START
    int foreWidth = GDALGetRasterXSize(foreData);
    int foreHeight = GDALGetRasterYSize(foreData);
    GDALRasterBand* foreBand = foreData->GetRasterBand(1);
    GDALDataType foreDataType=foreBand->GetRasterDataType();
    float *foreData32F=NULL;
    unsigned char *foreData8UC=NULL;
    
    if (foreDataType==GDT_Byte){
      foreData8UC = new unsigned char[foreHeight*foreWidth];
      foreBand->RasterIO(GF_Read, 0, 0,
			 foreWidth, foreHeight,
			 foreData8UC, foreWidth, foreHeight,
			 foreDataType, 0, 0);
    }
    if (foreDataType==GDT_Float32){
      foreData32F = new float[foreHeight*foreWidth];
      foreBand->RasterIO(GF_Read, 0, 0,
			 foreWidth, foreHeight,
			 foreData32F, foreWidth, foreHeight,
			 foreDataType, 0, 0);
    }
    //read the foreground - END
    
    //read the background - START
    int backWidth = GDALGetRasterXSize(backData);
    int backHeight = GDALGetRasterYSize(backData);
    GDALRasterBand* backBand = backData->GetRasterBand(1);
    GDALDataType backDataType=backBand->GetRasterDataType();
    
    float *backData32F=NULL;
    unsigned char *backData8UC=NULL;
    if (backDataType==GDT_Byte){
      backData8UC = new unsigned char[backHeight*backWidth];
      backBand->RasterIO(GF_Read, 0, 0,
			 backWidth, backHeight,
			 backData8UC, backWidth, backHeight,
			 backDataType, 0, 0);
    }
    if (backDataType==GDT_Float32){
      backData32F = new float[backHeight*backWidth];
      backBand->RasterIO(GF_Read, 0, 0,
			 backWidth, backHeight,
			 backData32F, backWidth, backHeight,
			 backDataType, 0, 0);
    }
    //read the background - END
    
    vector<double> forePixel;
    forePixel.resize(2);
    vector<double> backPixel;
    vector<double> geoCoord;
    float offset = 0;
    int numSamples = 0;
    
    for (int j = 0; j < foreHeight; j=j+40){
      for (int i = 0; i < foreWidth; i=i+40){
	
	forePixel[0] = i; forePixel[1] = j;
        
	geoCoord  = pixelToGeographic(forePixel, foreData);
	backPixel = geographicToPixel(geoCoord, backData);
        
	if ((backPixel[0]>0) && (backPixel[0]<backWidth) &&
	    (backPixel[1]>0) && (backPixel[1]<backHeight)){
	  
	  int x = (int)floor(backPixel[0]);
	  int y = (int)floor(backPixel[1]);
          
	  if ((foreDataType==GDT_Float32) &&
	      (backDataType==GDT_Byte) &&
	      !isnan(foreData32F[j*foreWidth+i]))
	    {
	      offset = offset + 255*foreData32F[j*foreWidth+i]-backData8UC[y*backWidth+x];
	      numSamples++;
	    }
	  if ((foreDataType == GDT_Float32) &&
	      (backDataType == GDT_Float32) &&
	      !isnan(foreData32F[j*foreWidth+i]))
	    {
	      //cout<<"assembler_util.cc:: computeOffset: test3.6"<<endl;
	      offset = offset + foreData32F[j*foreWidth+i]-backData32F[y*backWidth+x];
	      numSamples++;
	    }
	  if ((foreDataType == GDT_Byte) &&
	      (backDataType == GDT_Byte) &&
	      !isnan(foreData8UC[j*foreWidth+i]))
	    {
	      offset = offset + foreData8UC[j*foreWidth+i]-backData8UC[y*backWidth+x];
	      numSamples++;
	    }
	  
          
	} 
	geoCoord.resize(0);
	backPixel.resize(0);
        
      }
    }
    //cout<<"assembler_util.cc:: computeOffset: test4"<<endl;
    
    if (numSamples!=0){
      offset = offset/numSamples;
    }
    //cout<<"assembler_util.cc:: computeOffset: offset = "<<offset<<endl;
    
    forePixel.resize(0);
    if (backData32F!=NULL){
      delete backData32F;
    }
    if (foreData32F!=NULL){
      delete foreData32F;
    }
    if (backData8UC!=NULL){
      delete backData8UC;
    }
    if (foreData8UC!=NULL){
      delete foreData8UC;
    }
    
    return offset;
    
  }
  
  
  //return the set of filenames that overlap with data.
  vector<string>
  findOverlapFiles(vector<string> tileFilenameArray, GDALDataset* data)
  {
    double geoTransform[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    if (data->GetGeoTransform(geoTransform) == CE_Failure)
      cerr << "WARNING in selectTileFilename(): could not read geoTransform "
	   << "using:\n\t["
	   << geoTransform[0] << ", " << geoTransform[1] << ", "
	   << geoTransform[2] << ", " << geoTransform[3] << ", "
	   << geoTransform[4] << ", " << geoTransform[5] << "]" << endl; 
    int cols = GDALGetRasterXSize(data);
    int rows = GDALGetRasterYSize(data);
    float leftPoint  = geoTransform[0];
    float rightPoint = geoTransform[0] + geoTransform[1]*cols;
    float topPoint   = geoTransform[3];
    float lowPoint   = geoTransform[3] + geoTransform[5]*rows; 
    float centerX = 0.5*(rightPoint+leftPoint);
    float centerY = 0.5*(lowPoint+topPoint); 
    float width   = rightPoint-leftPoint;
    float height  = topPoint-lowPoint; 
  
    cout << "GeoUtils.cc:: selectTileFilename(): left = " << leftPoint
	 << ", right = " << rightPoint << ", top = " << topPoint << ", low = "
	 <<  lowPoint << endl;
    
    vector<string> selectTileFilenameArray;
    selectTileFilenameArray.resize(0);

    for (unsigned int i=0; i < tileFilenameArray.size(); i++){

      GDALDataset* tile = (GDALDataset *) GDALOpen(tileFilenameArray[i].c_str(), GA_ReadOnly);
      
      if (tile!=NULL){

	//determine the bounding box of the background tile - START
	double tileGeoTransform[6];
	tile->GetGeoTransform(tileGeoTransform);
            
	int tileCols = GDALGetRasterXSize(tile);
	int tileRows = GDALGetRasterYSize(tile);

	float tileLeftPoint  = tileGeoTransform[0];
	float tileRightPoint = tileGeoTransform[0] + tileGeoTransform[1]*tileCols;
	float tileTopPoint   = tileGeoTransform[3];
	float tileLowPoint   = tileGeoTransform[3] + tileGeoTransform[5]*tileRows;
	
	float tileCenterX = 0.5*(tileRightPoint+tileLeftPoint);
	float tileCenterY = 0.5*(tileLowPoint+tileTopPoint);
	float tileWidth   = tileRightPoint-tileLeftPoint;
	float tileHeight  = tileTopPoint-tileLowPoint;  
	//determine the bounding box of the background tile - BACK
      
	//if (fabs(a.x - b.x) * 2 < (a.width + b.width)) &&
	//   (fabs(a.y - b.y) * 2 < (a.height + b.height)){
      
	//condition to check the overlap between to bboxes
	if ( (2*fabs(tileCenterX-centerX) < tileWidth + width ) &&
	     (2*fabs(tileCenterY-centerY) < tileHeight+ height)) {
	  selectTileFilenameArray.push_back(tileFilenameArray[i]);
	  
	}
	GDALClose((GDALDatasetH)tile);
      }
    }
    
    cout<<"GeoUtils.cc:: findOverlapFiles(): num selected tiles="<<selectTileFilenameArray.size()<<endl;
    return selectTileFilenameArray;
  }
  
  //finds the overlapping bounding box between to GDALDatasets
  struct bbox findOverlapPixelBBox(GDALDataset* refDS, GDALDataset *currDS)
  {
    vector<double> refGeographicBBox = computeGeographicBBox(refDS);
    struct bbox refGeoBBox;
    refGeoBBox.xl=refGeographicBBox[0];
    refGeoBBox.yt=refGeographicBBox[3];//[2];
    refGeoBBox.xr=refGeographicBBox[1];
    refGeoBBox.yb=refGeographicBBox[2];//[3];
    //cout<<"refGeoBBox: xl="<<refGeoBBox.xl<<" ,yt="<<refGeoBBox.yt
    //	<<" ,xr="<<refGeoBBox.xr<<" ,yb="<<refGeoBBox.yb<<endl;
    struct bbox refPixelBBox = geographicToPixelBBox(refGeoBBox, refDS);
    
    //cout<<"refPixelBBox: xl="<<refPixelBBox.xl<<" ,yt="<<refPixelBBox.yt
    //	<<" ,xr="<<refPixelBBox.xr<<" ,yb="<<refPixelBBox.yb<<endl;
    vector<double> currGeographicBBox = computeGeographicBBox(currDS);
   
    struct bbox currGeoBBox;
    currGeoBBox.xl=currGeographicBBox[0];
    currGeoBBox.yt=currGeographicBBox[3];//[2];
    currGeoBBox.xr=currGeographicBBox[1];
    currGeoBBox.yb=currGeographicBBox[2];//[3];
    //cout<<"currGeoBBox: xl="<<currGeoBBox.xl<<" ,yt="<<currGeoBBox.yt
    //	<<" ,xr="<<currGeoBBox.xr<<" ,yb="<<currGeoBBox.yb<<endl;
    struct bbox currPixelBBox = geographicToPixelBBox(currGeoBBox, refDS);
    //cout<<"currPixelBBox: xl="<<currPixelBBox.xl<<" ,yt="<<currPixelBBox.yt
    //	<<" ,xr="<<currPixelBBox.xr<<" ,yb="<<currPixelBBox.yb<<endl;

    struct bbox overlapBBox = computeBoxIntersection(refPixelBBox, currPixelBBox);

    return overlapBBox;
  }
  
  vector<float> computeUnionGeographicBBox(vector<string> filenameArray) 
  {
    vector<float> unionBBoxGeographic;
    unionBBoxGeographic.resize(4);

    GDALDataset* ds = (GDALDataset *) GDALOpen(filenameArray[0].c_str(), GA_ReadOnly); 
    vector<double> bbox = computeGeographicBBox(ds);
    unionBBoxGeographic[0]=bbox[0];
    unionBBoxGeographic[1]=bbox[1];
    unionBBoxGeographic[2]=bbox[2];
    unionBBoxGeographic[3]=bbox[3];
    GDALClose(ds);
    
    for (unsigned int k = 0; k < filenameArray.size(); k++){
       GDALDataset* ds = (GDALDataset *) GDALOpen(filenameArray[k].c_str(), GA_ReadOnly); 
       vector<double> bbox = computeGeographicBBox(ds);
       cout<<"bbox"<<bbox[0]<<", "<<bbox[1]<<", "<<bbox[2]<<", "<<bbox[3]<<endl;
       if (bbox[0]<unionBBoxGeographic[0]){
	 unionBBoxGeographic[0]=bbox[0];
       }
       if (bbox[1]>unionBBoxGeographic[0]){
	 unionBBoxGeographic[1]=bbox[1];
       }
       if (bbox[2]<unionBBoxGeographic[2]){
	 unionBBoxGeographic[2]=bbox[2];
       }
       if (bbox[3]>unionBBoxGeographic[3]){
	 unionBBoxGeographic[3]=bbox[3];
       }
       //close GDALDatset
       GDALClose(ds);
    }
    cout<<"********************unionBBoxGeographic "
	<<unionBBoxGeographic[0]<<", "
	<<unionBBoxGeographic[1]<<", "
	<<unionBBoxGeographic[2]<<", "
	<<unionBBoxGeographic[3]<<endl;

    return unionBBoxGeographic;
  }


  
  vector<double>
  computeGeographicBBox(GDALDataset* data)
  {
    double geoTransform[6] = { 999.0, 999.0, 999.0, 999.0, 999.0, 999.0 };
    if (data->GetGeoTransform(geoTransform) == CE_Failure){
      cerr << "WARNING in geotof_utilcc:: computeGeographicBBox(): could not read geoTransform "
	   << "using:\n\t["
	   << geoTransform[0] << ", " << geoTransform[1] << ", "
	   << geoTransform[2] << ", " << geoTransform[3] << ", "
	   << geoTransform[4] << ", " << geoTransform[5] << "]" << endl;
    }

    int width = GDALGetRasterXSize(data);
    int height = GDALGetRasterYSize(data);
    double leftPoint = geoTransform[0];
    double rightPoint = geoTransform[0] + geoTransform[1]*width;
    double topPoint = geoTransform[3];
    double lowPoint = geoTransform[3] + geoTransform[5]*height;  

    /*
    cout << "geotif_util.cc:: computeGeographicBBox(): BBox projected is: "
	 <<"left = " << leftPoint << ", right = " << rightPoint << ", top = "
	 << topPoint << ", low = " << lowPoint << endl;
    */

    vector<double> lonlatBB;
    lonlatBB.resize(4);
    
    vector<double> projCoord;
    projCoord.resize(2);
    projCoord[0]=leftPoint;
    projCoord[1]=topPoint;
    
    vector<double> geoCoord;
    geoCoord = projectedToGeographic(projCoord, data);
    lonlatBB[0]=geoCoord[0];
    lonlatBB[3]=geoCoord[1];
    geoCoord.resize(0);
    
    projCoord[0]=rightPoint;
    projCoord[1]=lowPoint;
    geoCoord = projectedToGeographic(projCoord, data);
    lonlatBB[1]=geoCoord[0];
    lonlatBB[2]=geoCoord[1];

    //lon lefttop,  lon rightbottom, lat leftbottom, lat lefttop
    //minlon, maxlon, minlat, maxlat
    return lonlatBB;
  }

  struct bbox changePixelBBoxCoords(GDALDataset* src, GDALDataset* dst, struct bbox srcPixelBBox)
  {

    struct bbox projectedBBox = pixelToProjectedBBox(srcPixelBBox, src);    
    //determine the pixel bbox in the dst image
    struct bbox dstPixelBBox = projectedToPixelBBox(projectedBBox, dst);

    return dstPixelBBox;
  }

  struct bbox  projectedToPixelBBox(struct bbox projectedBBox, GDALDataset* data)
  {

    vector<double> topLeftProjCoord;
    topLeftProjCoord.resize(2);
    topLeftProjCoord[0]=projectedBBox.xl;
    topLeftProjCoord[1]=projectedBBox.yt;
  
    vector<double> bottomRightProjCoord;
    bottomRightProjCoord.resize(2);
    bottomRightProjCoord[0]=projectedBBox.xr;
    bottomRightProjCoord[1]=projectedBBox.yb;
        
    vector<double> topLeftPix=projectedToPixel(topLeftProjCoord, data);
    vector<double> bottomRightPix=projectedToPixel(bottomRightProjCoord, data);
 
    int width = GDALGetRasterXSize(data);
    int height =  GDALGetRasterYSize(data);
    
    //make sure foreROI is within the boundaries of the original foreData - START
    if (topLeftPix[0]<0){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning xl is negative"<<endl;
      topLeftPix[0]=0;
    }

    if (topLeftPix[0]>=width-1){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning xl outside the boundaries of the data"<<endl;
      topLeftPix[0]=width-1;
    }
    
    if (topLeftPix[1]<0){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning yt is negative"<<endl;
      topLeftPix[1]=0;
    }

    if (topLeftPix[0]>=height-1){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning yt outside the boundaries of the data"<<endl;
      topLeftPix[0]=height-1;
    }
    
    if (bottomRightPix[0]<0){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning xr is negative"<<endl;
      bottomRightPix[0]=0;
    }
    
    if (bottomRightPix[0]>=width-1){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning xr outside the boundaries of the data"<<endl;
      bottomRightPix[0]=width-1;
    }
    
    if (bottomRightPix[1]<0){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning yb is negative"<<endl;
      bottomRightPix[1]=0;
    }
    
    if (bottomRightPix[1]>=height-1){
      cout<<"geotif_util.cc:: projectedToPixelBBox(): Warning yb outside the boundaries of the data"<<endl;
      bottomRightPix[1]=height-1;
    }

    struct bbox pixelBBox;
    
    pixelBBox.xl=topLeftPix[0];
    pixelBBox.yt=topLeftPix[1];
    pixelBBox.xr=bottomRightPix[0];
    pixelBBox.yb=bottomRightPix[1];
    
    return pixelBBox;
  }
  
  struct bbox pixelToProjectedBBox(struct bbox pixelBBox, GDALDataset* data)
  {

    vector<double> topLeftPixCoord;
    topLeftPixCoord.resize(2);
    topLeftPixCoord[0]=pixelBBox.xl;
    topLeftPixCoord[1]=pixelBBox.yt;
	
    vector<double> bottomRightPixCoord;
    bottomRightPixCoord.resize(2);
    bottomRightPixCoord[0]=pixelBBox.xr;
    bottomRightPixCoord[1]=pixelBBox.yb;

    vector<double> topLeftProjCoord=pixelToProjected(topLeftPixCoord, data);
    /*
    cout<<"geotif_util.cc:: pixelToProjectedBBox(): topLeftProjCoord="<<topLeftProjCoord[0]
	<<", "<<topLeftProjCoord[1]<<endl;
    */
    vector<double> bottomRightProjCoord=pixelToProjected(bottomRightPixCoord, data);
    /*
    cout<<"geotif_util.cc:: pixelToProjectedBBox(): bottomRightProjCoord="<<bottomRightProjCoord[0]
	<<", "<<bottomRightProjCoord[1]<<endl;
    */
    struct bbox projectedBBox;
    projectedBBox.xl=topLeftProjCoord[0];
    projectedBBox.yt=topLeftProjCoord[1];
    projectedBBox.xr=bottomRightProjCoord[0];
    projectedBBox.yb=bottomRightProjCoord[1];
    
    return projectedBBox;
  }

  
  //determines the pixel bounding box for a specific lonlat bounding box
  struct bbox geographicToPixelBBox(struct bbox geographicBBox, GDALDataset* data)
  {
    vector<double> topLeftLonLat, bottomRightLonLat;
    topLeftLonLat.resize(2);
    topLeftLonLat[0]=geographicBBox.xl;
    topLeftLonLat[1]=geographicBBox.yt;
    
    bottomRightLonLat.resize(2);
    bottomRightLonLat[0]=geographicBBox.xr;
    bottomRightLonLat[1]=geographicBBox.yb;
        
    //determine the pixel in the background image corresponding to the topleft lonlat 
    vector<double> topLeftPix = geographicToPixel(topLeftLonLat, data);    
    vector<double> bottomRightPix = geographicToPixel(bottomRightLonLat, data);
    
    struct bbox pixelBB;
    pixelBB.xl = /*(float)*/topLeftPix[0];
    pixelBB.yt = /*(float)*/topLeftPix[1];
    pixelBB.xr = /*(float)*/bottomRightPix[0];
    pixelBB.yb = /*(float)*/bottomRightPix[1];
    
    return pixelBB;
  }
  
  //struct bbox pixelToGeographicBBox_(struct bbox pixelBBox, GDALDataset* data)
  struct bbox pixelToGeographicBBox(struct bbox pixelBBox, GDALDataset* data)
  {
 
    struct bbox geographicBBox;
    
    vector<double> topLeftPixel; topLeftPixel.resize(2);
    topLeftPixel[0]=pixelBBox.xl;
    topLeftPixel[1]=pixelBBox.yt;
    vector<double> bottomRightPixel; bottomRightPixel.resize(2);
    bottomRightPixel[0]=pixelBBox.xr;
    bottomRightPixel[1]=pixelBBox.yb;
 
    vector<double> topLeftPoint = pixelToGeographic(topLeftPixel, data);
    vector<double> bottomRightPoint = pixelToGeographic(bottomRightPixel, data);
    
    geographicBBox.xl=topLeftPoint[0];
    geographicBBox.yt=topLeftPoint[1];
    geographicBBox.xr=bottomRightPoint[0];
    geographicBBox.yb=bottomRightPoint[1];

    return geographicBBox;
  }

  //returns the mem buffer for the intersection between tile and an geo dataset
  //also upsamples the data.
  unsigned char* make8UCTile(GDALDataset* ds, float upsampleFactor, 
	                     struct TilingParams &tileParams, int bandIndex)
  {
    int resampledImgWidth  = GDALGetRasterXSize(ds)/upsampleFactor;
    int resampledImgHeight = GDALGetRasterYSize(ds)/upsampleFactor;
    GDALRasterBand* band   = ds->GetRasterBand(1);
    
    int tileImgWidth, tileImgHeight;
    tileImgWidth =tileParams.back_xr-tileParams.back_xl;
    tileImgHeight=tileParams.back_yb-tileParams.back_yt;
    
    unsigned char* ucharTileImg = new unsigned char[tileImgWidth*tileImgHeight];
    
    for (int i=0; i <tileImgWidth*tileImgHeight; i++){
      ucharTileImg[i] = 0;
    }
    //compute the intersection between the current tile and the entire image

    struct bbox tileBox;
    tileBox.xl=tileParams.back_xl;
    tileBox.yt=tileParams.back_yt;
    tileBox.xr=tileParams.back_xr;
    tileBox.yb=tileParams.back_yb;
    
    struct bbox imageBox;
    imageBox.xl=0;
    imageBox.yt=0;
    imageBox.xr=resampledImgWidth;
    imageBox.yb=resampledImgHeight;
    struct bbox intersectBox = computeBoxIntersection(tileBox, imageBox);
    
    if (isRectangle(intersectBox)){
      if ((intersectBox.xr-intersectBox.xl==tileImgWidth)&&(intersectBox.xr-intersectBox.yt==tileImgHeight)){
	//TODO:use the original tileParams values here.
	band->RasterIO(GF_Read,
		       intersectBox.xl*upsampleFactor, intersectBox.yt*upsampleFactor,
		       (intersectBox.xr-intersectBox.xl)*upsampleFactor,
		       (intersectBox.yb-intersectBox.yt)*upsampleFactor, 
		       ucharTileImg, tileImgWidth, tileImgHeight,
		       GDT_Byte, 0, 0 );
	
      }
      else{
	int tempImgWidth = intersectBox.xr-intersectBox.xl;
	int tempImgHeight = intersectBox.yb-intersectBox.yt;
	unsigned char *ucharTempImg = new unsigned char[tempImgWidth*tempImgHeight];
        
	//fill the tile with noData values 
	band->RasterIO(GF_Read,
		       intersectBox.xl*upsampleFactor, intersectBox.yt*upsampleFactor,
		       (intersectBox.xr-intersectBox.xl)*upsampleFactor,
		       (intersectBox.yb-intersectBox.yt)*upsampleFactor, 
		       ucharTempImg, tempImgWidth, tempImgHeight,
		       GDT_Byte, 0, 0 );
	
	for (int i = intersectBox.xl; i < intersectBox.xr; i++){
	  for (int j = intersectBox.yt; j < intersectBox.yb; j++){
	    int ii = i-tileParams.back_xl;
	    int jj = j-tileParams.back_yt; 
	    ucharTileImg[jj*tileImgWidth+ii]=ucharTempImg[(int)floor(j-intersectBox.yt)*tempImgWidth +
							  (int)floor((i-intersectBox.xl))];
	  }	   
	}

	delete ucharTempImg;
	
      }
    }
    
    return (ucharTileImg);
    
  }
  
  float* make32FTile(GDALDataset* ds, float upsampleFactor, 
		     struct TilingParams &tileParams, int bandIndex)
  {
    int resampledImgWidth  = GDALGetRasterXSize(ds)/upsampleFactor;
    int resampledImgHeight = GDALGetRasterYSize(ds)/upsampleFactor;
    GDALRasterBand* band = ds->GetRasterBand(1);
    
    int tileImgWidth, tileImgHeight;
    tileImgWidth =tileParams.back_xr-tileParams.back_xl;
    tileImgHeight=tileParams.back_yb-tileParams.back_yt;
     
    float* floatTileImg = new float[tileImgWidth*tileImgHeight];
    
    for (int i=0; i <tileImgWidth*tileImgHeight; i++){
      floatTileImg[i] = 0;
    }
    //compute the intersection between the current tile and the entire image

    struct bbox tileBox;
    tileBox.xl=tileParams.back_xl;
    tileBox.yt=tileParams.back_yt;
    tileBox.xr=tileParams.back_xr;
    tileBox.yb=tileParams.back_yb;
    
    struct bbox imageBox;
    imageBox.xl=0;
    imageBox.yt=0;
    imageBox.xr=resampledImgWidth;
    imageBox.yb=resampledImgHeight;
    struct bbox intersectBox = computeBoxIntersection(tileBox, imageBox);
    
    if (isRectangle(intersectBox)){
        if ((intersectBox.xr-intersectBox.xl==tileImgWidth)&&(intersectBox.yb-intersectBox.yt==tileImgHeight)){
	//TODO:use the original tileParams values here.
	band->RasterIO(GF_Read,
		       intersectBox.xl*upsampleFactor, intersectBox.yt*upsampleFactor,
		       (intersectBox.xr-intersectBox.xl)*upsampleFactor,
		       (intersectBox.yb-intersectBox.yt)*upsampleFactor, 
		       floatTileImg, tileImgWidth, tileImgHeight,
		       GDT_Float32, 0, 0 );
	
      }
      else{
	int tempImgWidth = intersectBox.xr-intersectBox.xl;
	int tempImgHeight = intersectBox.yb-intersectBox.yt;
	float *floatTempImg = new float[tempImgWidth*tempImgHeight];
        
	//fill the tile with noData values 
	band->RasterIO(GF_Read,
		       intersectBox.xl*upsampleFactor, intersectBox.yt*upsampleFactor,
		       (intersectBox.xr-intersectBox.xl)*upsampleFactor,
		       (intersectBox.yb-intersectBox.yt)*upsampleFactor, 
		       floatTempImg, tempImgWidth, tempImgHeight,
		       GDT_Float32, 0, 0 );
	
	for (int i = intersectBox.xl; i < intersectBox.xr; i++){
	  for (int j = intersectBox.yt; j < intersectBox.yb; j++){
	    int ii = i-tileParams.back_xl;
	    int jj = j-tileParams.back_yt; 
	    floatTileImg[jj*tileImgWidth+ii]=floatTempImg[(int)floor(j-intersectBox.yt)*tempImgWidth+
							  (int)floor((i-intersectBox.xl))];
	  }	   
	}

	delete floatTempImg;
	
      }
    }
    
    return (floatTileImg);
    
  }

  void saveGDALDataset(GDALDataset* data, float *data32F)
  {  
    int width = GDALGetRasterXSize(data);
    int height = GDALGetRasterYSize(data);
    GDALRasterBandH hBandout = GDALGetRasterBand(data, 1);
    CPLErr imError   = CE_None;
    imError = GDALRasterIO( hBandout, GF_Write, 0, 0, width, height,
			    data32F, width, height, GDT_Float32, 0, 0);
  }
  
  void saveGDALDataset(GDALDataset* data, unsigned char *data8UC)
  {
    int width = GDALGetRasterXSize(data);
    int height = GDALGetRasterYSize(data);
    GDALRasterBandH hBandout = GDALGetRasterBand(data, 1);
    CPLErr imError = CE_None;
    imError = GDALRasterIO( hBandout, GF_Write, 0, 0, width, height,
			    data8UC, width, height, GDT_Byte, 0, 0);
  }
  
  CPLErr saveGDALDataset(unsigned char *data8UC, int height, int width, int numBands,
			 double*geoTransform,  char *pszSRS_WKT, string filename)
  {
    
    GDALDriverH hDriver;
    GDALDatasetH hDstDS;
    char **pzOptions = NULL;
    
    cout << "GeoUtils.cc:: saveGDALDataset(): writing file " << filename << endl;
    
    GDALAllRegister();
    
    hDriver = GDALGetDriverByName("GTiff");
    if (hDriver==NULL){
    }
    
    pzOptions = CSLSetNameValue( pzOptions, "INTERLEAVE", "PIXEL" );
    
    hDstDS  = GDALCreate( hDriver, filename.c_str(),
			  width, height,
			  numBands, GDT_Byte , pzOptions );
    GDALSetGeoTransform( hDstDS, geoTransform );
      
    GDALSetProjection( hDstDS, pszSRS_WKT );
    
    CPLErr imError   = CE_None;
    GDALRasterBandH hBandout;
    if (numBands==1){ 
      hBandout= GDALGetRasterBand( hDstDS, 1 );
      imError = GDALRasterIO( hBandout, GF_Write, 0, 0, width, height,
			      data8UC, width, height, GDT_Byte, 0, 0);
    }
    if (numBands==3){ 
      unsigned char *dataBand = new unsigned char[width*height];
      
      for (int i=0; i < width*height; i++){
	dataBand[i] = data8UC[i];
      }
      hBandout= GDALGetRasterBand( hDstDS, 1 );
      imError = GDALRasterIO(hBandout, GF_Write, 0, 0, width, height,
			     dataBand, width, height, GDT_Byte, 0, 0);
      for (int i=0; i < width*height; i++){
	dataBand[i]=data8UC[i+width*height];
      }
      hBandout= GDALGetRasterBand( hDstDS, 2 );
      imError = GDALRasterIO(hBandout, GF_Write, 0, 0, width, height,
			     dataBand, width, height, GDT_Byte, 0, 0);
      
      for (int i=0; i < width*height; i++){
	dataBand[i]=data8UC[i+2*width*height];
      }
      hBandout= GDALGetRasterBand( hDstDS, 3 );
      imError = GDALRasterIO(hBandout, GF_Write, 0, 0, width, height,
			     dataBand, width, height, GDT_Byte, 0, 0);
      delete dataBand;
    }
    GDALClose( hDstDS );
    hBandout   = NULL;
    hDstDS   = NULL;
    
    return CE_None;
  }
  
  //TODO: allow this to write multiple channels
  CPLErr saveGDALDataset(float *data32F, int height, int width,
			 double*geoTransform, char *pszSRS_WKT, string filename)
  {
    cout << "GeoUtils.cc:: saveGDALDataset(): writing filename " << filename << endl;
    
    GDALDriverH hDriver;
    GDALDatasetH hDstDS;
    char **pzOptions = NULL;
    int imageBpp = 1;
    
    GDALAllRegister();
    
    hDriver = GDALGetDriverByName("GTiff");
    if (hDriver==NULL){cout<<"hDriver = NULL"<<endl;}
    
    pzOptions = CSLSetNameValue( pzOptions, "INTERLEAVE" , "PIXEL");
    //pzOptions = CSLSetNameValue( pzOptions, "COMPRESSION" , "LZW");
 
    int validOptions = GDALValidateCreationOptions (hDriver, pzOptions); 	
    
    hDstDS  = GDALCreate(hDriver, filename.c_str(),
			 width, height,
			 imageBpp, GDT_Float32, /*pzOptions*/NULL );

    GDALSetGeoTransform( hDstDS, geoTransform );

    GDALSetProjection( hDstDS, pszSRS_WKT );
    
    GDALRasterBandH hBandout = GDALGetRasterBand( hDstDS, 1 );
    CPLErr imError   = CE_None;
    
    imError = GDALRasterIO(hBandout, GF_Write, 0, 0, width, height,
			   data32F, width, height, GDT_Float32, 0, 0);
    
    GDALClose (hDstDS);
    hBandout = NULL;
    hDstDS   = NULL;
    
    return imError;
  } 
 
 
 
}
