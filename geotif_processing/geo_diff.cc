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

#include "../common/StringUtils.h"
#include "CoordTransform.h"

using namespace std;
using namespace at;

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) 
{

  std:: string filename1 = argv[1];
  std:: string filename2 = argv[2];
  std:: string resDirname = argv[3];
  
  GDALAllRegister();
     
  //read the fore data  
  GDALDataset* data1 = (GDALDataset *) GDALOpen(filename1.c_str(), GA_ReadOnly); 

  int width1  = GDALGetRasterXSize(data1);
  int height1 = GDALGetRasterYSize(data1);
  GDALRasterBand* bandData1 = data1->GetRasterBand(1);
  GDALDataType typeData1    =bandData1->GetRasterDataType();

  float * f32Data1 = NULL; 
  unsigned char *uC8Data1 = NULL;
  double noDataVal1;

  if (typeData1==GDT_Byte){ 
    uC8Data1 = new unsigned char[height1*width1];
    bandData1->RasterIO(GF_Read, 0, 0, 
		       width1, height1,  
		       uC8Data1, width1, height1, 
		       typeData1, 0, 0);
    int pbSuccess = 0;
    noDataVal1=0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData1, &pbSuccess);
    if (pbSuccess==1){
      noDataVal1= dSNoDataVal;
    } 	
  }

  if (typeData1==GDT_Float32){ 
    f32Data1 = new float[height1*width1];
    bandData1->RasterIO(GF_Read, 0, 0, 
		       width1, height1,  
		       f32Data1, width1, height1, 
		       typeData1, 0, 0);
    int pbSuccess = 0;
    noDataVal1=0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData1, &pbSuccess);
    if (pbSuccess==1){
      noDataVal1= dSNoDataVal;
    } 	
  }

  //read the background file
  GDALDataset* data2 = (GDALDataset *) GDALOpen(filename2.c_str(), GA_ReadOnly); 
  int width2  = GDALGetRasterXSize(data2);
  int height2 = GDALGetRasterYSize(data2);
  GDALRasterBand* bandData2 = data2->GetRasterBand(1);
  GDALDataType typeData2    =bandData2->GetRasterDataType();

  float * f32Data2 = NULL; 
  unsigned char *uC8Data2 = NULL;
  double noDataVal2;

  if (typeData2==GDT_Byte){ 
    uC8Data2 = new unsigned char[height2*width2];
    bandData2->RasterIO(GF_Read, 0, 0, 
		       width2, height2,  
		       uC8Data2, width2, height2, 
		       typeData2, 0, 0);
    int pbSuccess = 0;
    noDataVal2=0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData2, &pbSuccess);
    if (pbSuccess==2){
      noDataVal2= dSNoDataVal;
    } 	
  }

  if (typeData2==GDT_Float32){ 
    f32Data2 = new float[height2*width2];
    bandData2->RasterIO(GF_Read, 0, 0, 
		       width2, height2,  
		       f32Data2, width2, height2, 
		       typeData2, 0, 0);
    int pbSuccess = 0;
    noDataVal2=0;
    double dSNoDataVal=GDALGetRasterNoDataValue ((GDALRasterBandH)bandData2, &pbSuccess);
    if (pbSuccess==1){
      noDataVal2= dSNoDataVal;
    } 	
  }

  if (height1!=height2){
    std::cout<<"Error data height is different"<<endl;
    GDALClose(data1);
    GDALClose(data2);
    return 0;
  }
  if (width1!=width2){
    std::cout<<"Error data width is different"<<endl;
    GDALClose(data1);
    GDALClose(data2);
    return 0;
  }
  if (typeData1!=typeData2){
    std::cout<<"Error data types are different"<<endl;
    GDALClose(data1);
    GDALClose(data2);
    return 0;
  }
  
  float diff=0;
  float maxDiff = -1000000000.0;
  float minDiff =  1000000000.0;
  int numErrors=0;
  if ((typeData1==GDT_Float32)&&(typeData2 == GDT_Float32)){
    for (int j=0; j < height1; j++){
      cout<<"line "<<j<<" of "<<height1-1<<endl;
      for (int i=0; i < width1; i++){
	diff=f32Data1[i*width1+j]-f32Data2[i*width2+j];
        //cout<<"diff="<<diff<<endl;
        if (diff<minDiff){
	  minDiff=diff;
	}
	if (diff>maxDiff){
	  maxDiff=diff;
	}
        if (diff>1){
	  numErrors++;
	}
      }
    }
    cout<<"minDiff="<<minDiff<<", maxDiff="<<maxDiff<<", numErrors="<<(float)(numErrors*1.0)/(width1*height1)<<endl;
    GDALClose(data1);
    GDALClose(data2);
  }

  if ((typeData1 == GDT_Byte)&&(typeData2 == GDT_Byte)){
    for (int j=0; j < height1; j++){
      cout<<"line "<<j<<" of "<<height1-1<<endl;
      for (int i=0; i < width1; i++){
	diff=float(uC8Data1[i*width1+j]-uC8Data2[i*width2+j]);

        if (diff<minDiff){
	  minDiff=diff;
	}
	if (diff>maxDiff){
	  maxDiff=diff;
	}
        if (diff>1){
	  numErrors++;
	}
      }
    }
    cout<<"minDiff="<<minDiff<<", maxDiff="<<maxDiff<<", numErrors="<<(float)(numErrors*1.0)/(width1*height1)<<endl;
    GDALClose(data1);
    GDALClose(data2);
  }

  GDALClose(data1);
  GDALClose(data2);
 
  return 0;
}

















