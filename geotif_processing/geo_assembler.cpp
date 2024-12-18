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


#include <gdal_priv.h>

#include "../common/StringUtils.h"
#include "CoordTransform.h"
#include "GeoAssembler.h"
#include "GeoTiling.h"

using namespace std;
using namespace at;

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) 
{

  std::string foreFilename;
  std::string backFilename;
  std::string backDRGFilename;
  std::string foreDRGFilename;
  std::string resDir;
  std::string configFilename = "assembler_settings.txt";
  std::string mode;//it can be DEM or DRG


  int modeFound = 0;
  int backFilenameFound = 0;
  int foreFilenameFound = 0;
  int backDRGFilenameFound = 0;
  int foreDRGFilenameFound = 0;
  int resultsDirnameFound = 0;
  
  for (int i = 1; i < argc; i++) { //We will iterate over argv[] to get the parameters stored inside.
    //Note that we're starting on 1 because we don't need to know the
    //path of the program, which is stored in argv[0]
    if (i + 1 != argc) { // Check that we haven't finished parsing already

      if (string(argv[i]).compare("-m")==0) {
	mode = argv[i + 1];
	modeFound = 1;
      }
      else if (string(argv[i]).compare("-be")==0) {
	backFilename = argv[i + 1];
	backFilenameFound = 1;
      }
      else if (string(argv[i]).compare("-bp")==0) {
	backDRGFilename = argv[i + 1];
	backDRGFilenameFound = 1;
      }
      else if (string(argv[i]).compare("-fe")==0) {
	foreFilename = argv[i + 1];
	foreFilenameFound = 1;
      }
      else if (string(argv[i]).compare("-fp")==0) {
	foreDRGFilename = argv[i + 1];
	foreDRGFilenameFound = 1;
      }
      else if (string(argv[i]).compare("-r")==0) {
	resDir = argv[i+1];
	resultsDirnameFound = 1;
      }
      else if (string(argv[i]).compare("-c")==0) {
	configFilename = argv[i + 1];
      }    
    }
  }
 	
  int usageError = 0;
  
  if (modeFound == 0){
    cout<<"Error! Mode not found"<<endl;
    usageError = 1;
  }
  if (backFilenameFound==0) {
    cout<<"Error! Background filename not found"<<endl;
    usageError = 1;
  }
  if (backDRGFilenameFound==0) {
    cout<<"Error! Background Photo filename not found"<<endl;
    usageError = 1;
  }
  if (foreDRGFilenameFound==0) {
    cout<<"Error! Foreground Photo filename not found"<<endl;
    usageError = 1;
  }
  if (resultsDirnameFound==0) {
    cout<<"Error! result directory name not found"<<endl;
    usageError = 1;
  }
  
  if (usageError == 1){
    std::cout << "Usage is "<< argv[0]<< " -m mode(Elevation or Photo) -be backElevationFilename -bp backPhotoFilename"
                 "-fe foreElevationFilename -fp forePhotoFilename -r resultsDirname [-c configurationFilename]"<<endl;
    // Inform the user of how to use the program
    std::cin.get();
    exit(0);
  }
  
    
  cout.precision(12);
  struct AssemblerParams assemblerParams;
  ReadAssemblerConfigFile(configFilename, &assemblerParams);
  PrintSettings(&assemblerParams);
  
  GDALAllRegister();
 
  if ((mode.compare("DEM")==0)){
    //this should be set from the config or command line
    int overwriteFlag = 1;
    std::vector<vector<string> > assembledDEMTileFilenameArray = assembleGeoData(foreFilename, backFilename,
										 assemblerParams, "DEM", resDir, overwriteFlag);
  }
  
  if ((mode.compare("DRG")==0)){
    //this should be set from the config or command line
    int overwriteFlag = 1;
    std::vector<vector<string> >assembledDRGTileFilenameArray = assembleGeoData(foreDRGFilename, backDRGFilename,
										assemblerParams, "DRG", resDir, overwriteFlag);
  }
  
  if (mode.compare("DEM_DRG")==0){
    int overwriteFlag = 1;
    vector<struct TerrainFilename>assembledTileFilenameArray =  assembleGeoData(foreFilename, backFilename,
										foreDRGFilename, backDRGFilename,
										assemblerParams, resDir, overwriteFlag);
  }

  
  return 0;
}

















