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
#include "GeoTiling.h"


using namespace std;
using namespace at;

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) 
{

    string inputFilename;
    string outputDirname;
    string configFilename;

    int horTile, verTile, numPyrLevels;
    
    int inputFilenameFound = 0;
    int outputDirnameFound = 0;
    int horizontalTileFound = 0;
    int verticalTileFound = 0;
    int numPyramidLevelsFound = 0;
    
    for (int i = 1; i < argc; i++) { //We will iterate over argv[] to get the parameters stored inside.
      //Note that we're starting on 1 because we don't need to know the
      //path of the program, which is stored in argv[0]
      if (i + 1 != argc) { // Check that we haven't finished parsing already
	
	if (string(argv[i]).compare("-i")==0) {
	  inputFilename = argv[i + 1];
	  inputFilenameFound = 1;
	}
	else if (string(argv[i]).compare("-o")==0) {
	  outputDirname = argv[i + 1];
	  outputDirnameFound = 1;
	}
	else if (string(argv[i]).compare("-x")==0) {
	  horTile = atoi(argv[i + 1]);
	}
	else if (string(argv[i]).compare("-y")==0) {
	  verTile = atoi(argv[i + 1]);
	}
	else if (string(argv[i]).compare("-n")==0) {
	  numPyrLevels = atoi(argv[i + 1]);
	}
	else if (string(argv[i]).compare("-c")==0) {
	  configFilename = argv[i + 1];
	}    
      }
    }
    
  	
  int usageError = 0;
  
  if (inputFilenameFound == 0){
    cout<<"Error! Input filename not found"<<endl;
    usageError = 1;
  }
  
  if (outputDirnameFound==0) {
    cout<<"Error! Output dirname not found"<<endl;
    usageError = 1;
  }

  if (usageError == 1){
    std::cout << "Usage is "<< argv[0]<< " -i inputFilename -o outputDirname "
              <<"-x horizontalTile -y verticalTile -n numPyramidLevels [-c configurationFilename]"<<endl;
    // Inform the user of how to use the program
    std::cin.get();
    exit(0);
  }
    
    //resampling factor: < 1 is downsampling the image, > 1 is upsampling the image   
    
    struct ResampleParams resampleParams;
    int error = ReadResampleConfigFile(configFilename, &resampleParams);
    PrintResampleParams(&resampleParams);
     
   cout<<"pyrtiles_geotiff: outputDirname="<<outputDirname<<endl;

   if (resampleParams.imageType == 0){//DEM
     
     GDALAllRegister();
     // Assume a geo tiff in easting northing
     GDALDataset* dataset = (GDALDataset *) GDALOpen(inputFilename.c_str(), GA_ReadOnly); 
     
     //create the outputdirname
     string mkdirCommand;
     mkdirCommand="mkdir "+outputDirname;
     system(mkdirCommand.c_str());
     cout<<mkdirCommand<<endl;

     string elevationDirname=outputDirname+"/Elevation";
     mkdirCommand.clear();
     mkdirCommand="mkdir "+elevationDirname;
     cout<<mkdirCommand<<endl;
     system(mkdirCommand.c_str());

     std::vector<std::vector<struct TilingParams> > demPyrTileParams;

     vector<string> filenameArray;
     filenameArray.push_back(inputFilename);
     vector<float> unionBBoxGeographic = computeUnionGeographicBBox(filenameArray);
     struct bbox unionGeoBBox;
     unionGeoBBox.xl = unionBBoxGeographic[0];
     unionGeoBBox.yb = unionBBoxGeographic[2];
     unionGeoBBox.xr = unionBBoxGeographic[1];
     unionGeoBBox.yt = unionBBoxGeographic[3];
     
     demPyrTileParams = makeVirtualTiles(dataset, unionGeoBBox, resampleParams);

     writePyramidTiles(dataset, 0, 1, demPyrTileParams, "DEM", elevationDirname) ;
   }
   
   if (resampleParams.imageType == 1){//DRG 
 
      GDALAllRegister();
      // Assume a geo tiff in easting northing
      GDALDataset* dataset = (GDALDataset *) GDALOpen(inputFilename.c_str(), GA_ReadOnly); 
      
      //create the outputdirname
      string mkdirCommand;
      mkdirCommand="mkdir "+outputDirname;
      cout<<mkdirCommand<<endl;
      system(mkdirCommand.c_str());
      
      string photoDirname=outputDirname+"/Photo";
      mkdirCommand.clear();
      mkdirCommand="mkdir "+photoDirname;
      cout<<mkdirCommand<<endl;
      system(mkdirCommand.c_str());

      std::vector<std::vector<struct TilingParams> > drgPyrTileParams;

      vector<string> filenameArray;
      filenameArray.push_back(inputFilename);
      vector<float> unionBBoxGeographic = computeUnionGeographicBBox(filenameArray);
      struct bbox unionGeoBBox;
      unionGeoBBox.xl = unionBBoxGeographic[0];
      unionGeoBBox.yb = unionBBoxGeographic[2];
      unionGeoBBox.xr = unionBBoxGeographic[1];
      unionGeoBBox.yt = unionBBoxGeographic[3];
      drgPyrTileParams = makeVirtualTiles(dataset, unionGeoBBox, resampleParams);

      writePyramidTiles(dataset, 0, 1, drgPyrTileParams, "DRG", photoDirname) ;
   } 


   return 0;
}

















