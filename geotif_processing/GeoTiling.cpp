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
#include <ogr_spatialref.h>

#include "GeoTiling.h"

using namespace std;

namespace at 
{

 
  int ReadResampleConfigFile(string resampleConfigFilename, struct ResampleParams *resampleParams)
  {
    ifstream configFile (resampleConfigFilename.c_str());
    std::string line;
    double val; 
    std::string identifier;
    
    if (configFile.is_open()){ 
      
      cout << "Reading settings file "<< resampleConfigFilename<<"."<<endl; 
      
      while (!configFile.eof()){
	std::getline(configFile, line);
	
	if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines
	  
	  stringstream lineStream(line);
	  
	  string keyword("");
	  lineStream >> keyword;
	  
	  if (keyword.compare(string("IMAGE_TYPE"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->imageType = val;
	  }
	  if (keyword.compare(string("NUM_PYRAMID_LEVELS"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->numPyrLevels = val;
	  }
	  if (keyword.compare(string("PYRAMID_RESAMPLE_FACTOR"))==0){
	    float val; 
	    lineStream >> val;
	    resampleParams->pyrResampleFactor = val;
	  }
	  
	  if (keyword.compare(string("IMAGE_ROI_LEFT"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->imageROILeft=val;
	  }
	  if (keyword.compare(string("IMAGE_ROI_TOP"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->imageROITop=val;
	  }
	  if (keyword.compare(string("IMAGE_ROI_RIGHT"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->imageROIRight=val;
	  }
	  if (keyword.compare(string("IMAGE_ROI_BOTTOM"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->imageROIBottom=val;
	  }
	  
	  if (keyword.compare(string("TILE_SIZE"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->tileSize=val;
	  }
	  
	  if (keyword.compare(string("TILE_PADDING_TOP"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->tilePaddingTop=val;
	  }
	  
	  if (keyword.compare(string("TILE_PADDING_LEFT"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->tilePaddingLeft=val;
	  }
	  
	  if (keyword.compare(string("TILE_PADDING_RIGHT"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->tilePaddingRight=val;
	  }
	  
	  if (keyword.compare(string("TILE_PADDING_BOTTOM"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->tilePaddingBottom=val;
	  }
	  
	  if (keyword.compare(string("PYRAMID_MODE"))==0){
	    int val; 
	    lineStream >> val;
	    resampleParams->pyramidMode=val;
	  }
	  
	}
	
      }
    }
    
    return 0;
  }
  
  
  void PrintResampleParams(struct ResampleParams *resampleParams)
  {
    cout<<"IMAGE_TYPE = "<<resampleParams->imageType<<endl;
    cout<<"NUM_PYRAMID_LEVELS = "<<resampleParams->numPyrLevels<<endl;
    cout<<"PYRAMID_RESAMPLE_FACTOR = "<<resampleParams->pyrResampleFactor<<endl;
    cout<<"IMAGE_ROI_LEFT = "<<resampleParams->imageROILeft<<endl; 
    cout<<"IMAGE_ROI_TOP = "<<resampleParams->imageROITop<<endl;
    cout<<"IMAGE_ROI_RIGHT = "<<resampleParams->imageROIRight<<endl;
    cout<<"IMAGE_ROI_BOTTOM = "<<resampleParams->imageROIBottom<<endl;
    cout<<"TILE_SIZE = "<<resampleParams->tileSize<<endl;
    cout<<"TILE_PADDING_TOP = "<<resampleParams->tilePaddingTop<<endl;
    cout<<"TILE_PADDING_LEFT = "<<resampleParams->tilePaddingLeft<<endl;
    cout<<"TILE_PADDING_RIGHT = "<<resampleParams->tilePaddingRight<<endl;
    cout<<"TILE_PADDING_BOTTOM = "<<resampleParams->tilePaddingBottom<<endl;
    cout<<"PYRAMID_MODE = "<<resampleParams->pyramidMode<<endl;
  }

  std::vector<std::vector<struct TilingParams> > makeVirtualTiles(GDALDataset *data,
								  struct bbox unionGeoBBox,
								  ResampleParams resampleParams)
  {     
    vector<float> tilePaddingParams;
    tilePaddingParams.resize(4);
    tilePaddingParams[0]=resampleParams.tilePaddingTop;
    tilePaddingParams[1]=resampleParams.tilePaddingLeft;
    tilePaddingParams[2]=resampleParams.tilePaddingBottom;
    tilePaddingParams[3]=resampleParams.tilePaddingRight;
    
    struct TileConfig tileConfigParams;
    tileConfigParams.imageROILeft = resampleParams.imageROILeft;
    tileConfigParams.imageROITop = resampleParams.imageROITop;
    tileConfigParams.imageROIRight = resampleParams.imageROIRight;
    tileConfigParams.imageROIBottom = resampleParams.imageROIBottom;
    tileConfigParams.tileWidth = resampleParams.tileSize;
    tileConfigParams.tileHeight = resampleParams.tileSize;
    tileConfigParams.overlapHor = 0;
    tileConfigParams.overlapVer = 0;
    tileConfigParams.paddingTop = resampleParams.tilePaddingTop;
    tileConfigParams.paddingLeft = resampleParams.tilePaddingLeft;
    tileConfigParams.paddingRight = resampleParams.tilePaddingBottom;
    tileConfigParams.paddingBottom = resampleParams.tilePaddingRight;
    tileConfigParams.pyrResampleFactor = resampleParams.pyrResampleFactor;

    /*
    struct vector<float> unionGeographicBBox;
    unionGeographicBBox.resize(4);
    unionGeographicBBox[0]=unionGeoBBox.xl;
    unionGeographicBBox[1]=unionGeoBBox.xr;
    unionGeographicBBox[2]=unionGeoBBox.yb;
    unionGeographicBBox[3]=unionGeoBBox.yt;
    struct bbox unionBBoxPixel = geographicToPixelBBox(unionGeographicBBox, data);
    */
    struct bbox unionBBoxPixel = geographicToPixelBBox(unionGeoBBox, data);
    
    /*
    unionBBoxPixel.xl = round(unionBBoxPixel.xl/16)*16;
    unionBBoxPixel.yt = round(unionBBoxPixel.yt/16)*16;
    unionBBoxPixel.xr = round(unionBBoxPixel.xr/16)*16;
    unionBBoxPixel.yb = round(unionBBoxPixel.yb/16)*16;
    */
    //TODO: round the unionBBoxPixel to integer numbers.
    cout<<"pyrtiles_geotif_util.cc:: MakeTileParams(): unionBBoxPixel"
	<<": xl="<<unionBBoxPixel.xl
    	<<", yt="<<unionBBoxPixel.yt
      	<<", xr="<<unionBBoxPixel.xr
      	<<", yb="<<unionBBoxPixel.yb<<endl;

    std::vector<std::vector<struct TilingParams> >pyrTileParamsArray;

    Tiling tiler = Tiling(tileConfigParams);
    tiler.printTilingParams();
    tiler.processQuadTree(unionBBoxPixel,resampleParams.pyramidMode);
    pyrTileParamsArray = tiler.m_pyrTileParamsArray;
    //exit(1);
    
    /*
    vector<double> pixel; pixel.resize(2);    
    pixel[0]=pyrTileParamsArray[0][0].back_xl;
    pixel[1]=pyrTileParamsArray[0][0].back_yt;
    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): pixel_topleft="<<pixel[0]<<", "<<pixel[1]<<endl;
    double adfGeoTransform[6];
    data->GetGeoTransform(adfGeoTransform);
    vector<double> projected = pixelToProjected(pixel, data);
    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): DATA RESOLUTION="<<adfGeoTransform[1]<<", "<<adfGeoTransform[5]<<endl;

    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): Padded Terrain  TIE POINT="<<projected[0]<<", "<<projected[1]<<endl;

    pixel[0]=0; pixel[1]=0;
    projected = pixelToProjected(pixel, data);
    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): HIRISE TIE POINT="<<projected[0]<<", "<<projected[1]<<endl;

    //HIRISE
    projected[0]= 8131125.0;
    projected[1]=-259188.0;
    pixel = projectedToPixel(projected, data);
    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): HIRISE topleft pixel in padded terrain="
	<<pixel[0]<<", "<<pixel[1]<<endl;

    //HRSC
    projected[0]= 8066131.994;
    projected[1]=-200891.239;
    pixel = projectedToPixel(projected, data);
    cout<<"pyrTiles_geotif_util.cc:: makeVirtualTiles(): HRSC topleft pixel in padded terrain="
	<<pixel[0]<<", "<<pixel[1]<<endl;

    //exit(1);
    */
    return pyrTileParamsArray;
    
  }

  void checkVirtualTiles(std::vector<std::vector<struct TilingParams> >pyrTileArray, int tileSize)
  {
    cout<<"pyrtiles_geotif_util.cc::checkVirtualTiles(): numPyramidLevels="<<pyrTileArray.size()<<endl;

    int numWrongTiles=0;

    for (int l=0; l<pyrTileArray.size(); l++){
      for (int t=0; t<pyrTileArray[l].size(); t++){
	if ((pyrTileArray[l][t].back_xr-pyrTileArray[l][t].back_xl)!=tileSize){
	  cout<<"pyrtiles_geotif_util.cc::checkVirtualTiles():Incorrect tile at"
	      <<" l="<<l
	      <<", horIndex="<<pyrTileArray[l][t].horTileIndex
	      <<", verIndex="<<pyrTileArray[l][t].verTileIndex
	      <<", "<<pyrTileArray[l][t].back_xr-pyrTileArray[l][t].back_xl<<", "<<tileSize<<endl;
	  numWrongTiles++;
	}	
      }
    }
    if (numWrongTiles>0){
      cout<<"numWrongTiles="<<numWrongTiles<<endl;
      exit(1);
    }
  }
  
  //computes a new TilingParams strategy (for drg) using an existing TilingParams strategy (for dem)
  std::vector<std::vector<struct TilingParams> > changeVirtualTiles(GDALDataset* srcDS, struct ResampleParams srcResampleParams,
								    GDALDataset* dstDS, struct ResampleParams dstResampleParams,
								    std::vector<std::vector<struct TilingParams> >srcPyrTileParamsArray)
  {
    
    std::vector<double> topLeftPixelDEM;
    topLeftPixelDEM.resize(2);
    std::vector<double> bottomRightPixelDEM;
    bottomRightPixelDEM.resize(2);
    
    std::vector<std::vector<struct TilingParams> >dstPyrTileParamsArray;
    dstPyrTileParamsArray.resize(srcPyrTileParamsArray.size());
    
    float resampleFactor = 1;
    
    for (int p = 0; p <srcPyrTileParamsArray.size(); p++){
      
      dstPyrTileParamsArray[p].resize(srcPyrTileParamsArray[p].size());
      
      for (int m = 0; m <srcPyrTileParamsArray[p].size(); m++){
	
	topLeftPixelDEM[0] = srcPyrTileParamsArray[p][m].back_xl/resampleFactor;
	topLeftPixelDEM[1] = srcPyrTileParamsArray[p][m].back_yt/resampleFactor;

	//compute the projected coordinates in full res DEM DS
	std::vector<double> topLeftProjected = pixelToProjected(topLeftPixelDEM, srcDS);
	//compute the pixel coordinates in the full res DRG DS
	std::vector<double> topLeftPixelDRG  = projectedToPixel(topLeftProjected, dstDS);
	//append the DRG tile padding
	topLeftPixelDRG[0] = topLeftPixelDRG[0]*resampleFactor-dstResampleParams.tilePaddingLeft;   
	topLeftPixelDRG[1] = topLeftPixelDRG[1]*resampleFactor-dstResampleParams.tilePaddingTop;
	
	bottomRightPixelDEM[0] = (srcPyrTileParamsArray[p][m].back_xr-srcResampleParams.tilePaddingRight)/resampleFactor;
	bottomRightPixelDEM[1] = (srcPyrTileParamsArray[p][m].back_yb-srcResampleParams.tilePaddingBottom)/resampleFactor;
	std::vector<double> bottomRightProjected;
	bottomRightProjected = pixelToProjected(bottomRightPixelDEM, srcDS);
	std::vector<double> bottomRightPixelDRG;
	bottomRightPixelDRG = projectedToPixel(bottomRightProjected, dstDS);
	//append the DRG tile padding
	bottomRightPixelDRG[0] = bottomRightPixelDRG[0]*resampleFactor+dstResampleParams.tilePaddingRight;
	bottomRightPixelDRG[1] = bottomRightPixelDRG[1]*resampleFactor+dstResampleParams.tilePaddingBottom;

	/*
	bottomRightPixelDRG.resize(2);
        bottomRightPixelDRG[0]=topLeftPixelDRG[0]+drgResampleParams.tilePaddingLeft+drgResampleParams.tilePaddingRight+drgResampleParams.tileSize;
	bottomRightPixelDRG[1]=topLeftPixelDRG[1]+drgResampleParams.tilePaddingTop+drgResampleParams.tilePaddingBottom+drgResampleParams.tileSize;
	*/ 

	dstPyrTileParamsArray[p][m].back_xl = round(topLeftPixelDRG[0]);
	dstPyrTileParamsArray[p][m].back_yt = round(topLeftPixelDRG[1]);
	dstPyrTileParamsArray[p][m].back_xr = round(bottomRightPixelDRG[0]);
	dstPyrTileParamsArray[p][m].back_yb = round(bottomRightPixelDRG[1]);
	/*
	  if (p>=6){
	  cout<<"p="<<p<<", m="
	  <<m<<", xl="<<outPyrTileParamsArray[p][m].back_xl<<", yt="<<outPyrTileParamsArray[p][m].back_yt
	  <<", back_xl="<<topLeftPixelDRG[0]<<", back_yt="<<topLeftPixelDRG[1]<<endl;
	  }
	*/
	dstPyrTileParamsArray[p][m].horTileIndex=srcPyrTileParamsArray[p][m].horTileIndex;
	dstPyrTileParamsArray[p][m].verTileIndex=srcPyrTileParamsArray[p][m].verTileIndex;
	/*
	  if (m==0){
	  cout<<"p="<<p<<", m="<<m<<",xl="<<outPyrTileParamsArray[p][m].back_xl
	  <<", yt="<<outPyrTileParamsArray[p][m].back_yt
	  <<", xr="<<outPyrTileParamsArray[p][m].back_xr
	  <<", yb="<<outPyrTileParamsArray[p][m].back_yb<<endl;
	  cout<<"p="<<p<<", m="<<m<<",xl="<<pyrTileParamsArray[p][m].back_xl
	  <<", yt="<<pyrTileParamsArray[p][m].back_yt
	  <<", xr="<<pyrTileParamsArray[p][m].back_xr
	  <<", yb="<<pyrTileParamsArray[p][m].back_yb<<endl;
	  cout<<"----------------------------------------"<<endl;
	  }
	*/
      }
      
      resampleFactor = resampleFactor*0.5;
    }
    
    return dstPyrTileParamsArray;
  }


  
  void writePyramidTiles(GDALDataset* origDataset, int bandIndex, int numBands, 
		         std::vector<std::vector<struct TilingParams> >&pyrTileParamsArray,
		         string mode, string outputDirname) 
  {
    
    int numPyrLevels = pyrTileParamsArray.size()-1;
    
    double resampleFactor = 1.0;
    

    int tileImgWidth, tileImgHeight;


    //allocate memory for each tile - END
    
    for (int p = 0; p <= numPyrLevels; p++){
      
      cout<<"pyrLevel="<<p<<endl;
         
      int k,l;

      cout<<"numTiles = "<<pyrTileParamsArray[p].size()<<endl;
      //downsample the data
      GDALDatasetH hDownsample = resampleGDALDataset(origDataset, resampleFactor,
      						     "downsample.tif", string("round"), 0);
      
      for (unsigned int m = 0; m <pyrTileParamsArray[p].size(); m++){

	
	stringstream horTileIndex;
	horTileIndex<<pyrTileParamsArray[p][m].horTileIndex;
	stringstream verTileIndex;
	verTileIndex<<pyrTileParamsArray[p][m].verTileIndex;
	stringstream pyrLevel;
	pyrLevel<<p;

	//if ((p==0)&&
	//    (pyrTileParamsArray[p][m].horTileIndex==139)&&
	//    (pyrTileParamsArray[p][m].verTileIndex==92))
	//  {
	struct bbox tilePixelBBox;
	tilePixelBBox.xl=pyrTileParamsArray[p][m].back_xl;
	tilePixelBBox.yt=pyrTileParamsArray[p][m].back_yt;
	tilePixelBBox.xr=pyrTileParamsArray[p][m].back_xr;
	tilePixelBBox.yb=pyrTileParamsArray[p][m].back_yb;

	if (mode.compare("DEM")==0){
	   
	  float* floatTileImg = make32FTile(origDataset, (float)(1.0/resampleFactor), 
					    pyrTileParamsArray[p][m], bandIndex);

	  string tileFilename =  outputDirname+"/dem_"+pyrLevel.str()+"_"+horTileIndex.str()+"_"+verTileIndex.str()+".tif";
	  cout<<"tileFilename="<<tileFilename<<endl;
	  GDALDatasetH hTileDS = resizeGDALDataset((GDALDataset*)hDownsample, tilePixelBBox, tileFilename);
	  CPLErr err = GDALReprojectImage (origDataset, origDataset->GetProjectionRef(), hTileDS, NULL,
	  				   GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL); 
	  
	  saveGDALDataset((GDALDataset*)hTileDS, floatTileImg);
	  
	  GDALClose(hTileDS);
	  delete floatTileImg;

	}
	
	if (mode.compare("DRG")==0){
	  
	  unsigned char* ucharTileImg = make8UCTile(origDataset, (float)(1.0/resampleFactor), 
						    pyrTileParamsArray[p][m], bandIndex);
	  
	  string tileFilename =  outputDirname+"/drg_"+pyrLevel.str()+"_"+horTileIndex.str()+"_"+verTileIndex.str()+".tif";
	  cout<<"tileFilename="<<tileFilename<<endl;
	  GDALDatasetH hTileDS = resizeGDALDataset((GDALDataset*)hDownsample, tilePixelBBox, tileFilename);
	  CPLErr err = GDALReprojectImage (origDataset, origDataset->GetProjectionRef(), hTileDS, NULL,
	  				   GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL);
	  
	  saveGDALDataset((GDALDataset*)hTileDS, ucharTileImg);
	  GDALClose( hTileDS);
	  delete ucharTileImg;

	}
	//	  }
      }//m - tile index
      
      resampleFactor = resampleFactor*0.5;//resampleParams.pyrResampleFactor;
      GDALClose(hDownsample);
    }//p - pyr level

  }  

}
