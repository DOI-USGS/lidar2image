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

#include "../common/StringUtils.h"
#include "../common/FileListUtils.h"
#include "../common/ImageProcessing.h"

#include "../geotif_processing/CoordTransform.h"
#include "../geotif_processing/GeoAssembler.h"
#include "../geotif_processing/GeoUtils.h"
#include "../geotif_processing/GeoTiling.h"


using namespace std;

namespace at 
{

  //creates Antares pyramids of assembled tiles
  //returns the filenames of the assembled tiles.
  vector<struct TerrainFilename> pyrAssembleGeoData(GDALDataset *foreDEM, GDALDataset *backDEM,
						    GDALDataset *foreDRG, GDALDataset *backDRG,
						    struct bbox unionGeoBBox,
						    float offsetDEM, struct AssemblerParams assemblerParams, 
						    string thisResDir, int overwriteFlag)
  {
    vector<struct TerrainFilename> assembledTerrainFilenameArray;
    vector<vector<string> >assembledDEMFilenameArray;
    vector<vector<string> >assembledDRGFilenameArray;
    
    ResampleParams demResampleParams;
    demResampleParams.imageType=0;
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
    drgResampleParams.imageType=1;
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

    
    cout<<"pyramid_assembler_antares::pyrAssembleGeoDataAntares(): Making Virtual Reference Layers..."<<endl;
    float foreToReferenceResampleFactorDEM;
    GDALDatasetH hReferenceDEM = makeReferenceLayer(foreDEM, backDEM,
                                                    unionGeoBBox, assemblerParams, "DEM",
                                                    &foreToReferenceResampleFactorDEM);
    float foreToReferenceResampleFactorDRG;
    GDALDatasetH hReferenceDRG = makeReferenceLayer(foreDRG, backDRG,
                                                    unionGeoBBox, assemblerParams, "DRG", 
						    &foreToReferenceResampleFactorDRG);
    //multi-scale tiling DEM and DRG strategy for the reference layer - END
    
    cout<<"pyramid_assembler_antares::pyrAssembleGeoDataAntares(): Making Virtual Tiles..."<<endl;
    //multi-scale tiling DEM and DRG strategy for the reference layer - START
    //TODO: check that coverage for unionBBox and hReferenceDEM is the same.
    //unionBBoxGeographic is the coords of the back and hReferenceDEM in the coordinates of the fore.
    std::vector<std::vector<struct TilingParams> > pyrTileStrategyRefDEM;
    pyrTileStrategyRefDEM = makeVirtualTiles((GDALDataset*)hReferenceDEM,
					     unionGeoBBox,
					     demResampleParams);


    checkVirtualTiles(pyrTileStrategyRefDEM, demResampleParams.tileSize+1);

    std::vector<std::vector<struct TilingParams> > pyrTileStrategyRefDRG;
    pyrTileStrategyRefDRG = changeVirtualTiles((GDALDataset*)hReferenceDEM, demResampleParams,
                                               (GDALDataset*)hReferenceDRG, drgResampleParams,
					       pyrTileStrategyRefDEM);

    checkVirtualTiles(pyrTileStrategyRefDRG, drgResampleParams.tileSize+3);

    
    cout<<"pyramid_assembler_antares::pyrAssembleGeoDataAntares(): Making Data Reference Layers..."<<endl;
    //create the reference DEM data layer
    float foreToReferenceOffsetDEM;
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

    //create the reference DRG data layer
    float foreToReferenceOffsetDRG = 0;
    CPLErr err_drg = GDALReprojectImage (backDRG, backDRG->GetProjectionRef(), (GDALDataset*)hReferenceDRG, NULL,
      				         GRA_Bilinear, 0.0, 0.1,  NULL, NULL, NULL);

   
    cout<<"Making the Antares Pyramid ..."<<endl;
    //DEM and DRG background pyramid - START

    //DEM background pyramid    
    writePyramidTiles( (GDALDataset*)hReferenceDEM, 0,1, 
		       pyrTileStrategyRefDEM,
		       "DEM", thisResDir+"/Pyramid");
    
    //DRG background pyramid
    writePyramidTiles( (GDALDataset*)hReferenceDRG, 0,1, 
		       pyrTileStrategyRefDRG,
		       "DRG", thisResDir+"/Pyramid");
    
    //DEM and DRG  background pyramid - END
    
    cout<<"Making the assembled tiles ..."<<endl;
    //assembler - START
    
    assembledDEMFilenameArray = writeAssembledTiles(foreDEM, backDEM, hReferenceDEM,
						    foreToReferenceOffsetDEM,  assemblerParams.weightingMode, 
						    foreToReferenceResampleFactorDEM, overwriteFlag,
						    assemblerParams.backReference,
						    pyrTileStrategyRefDEM, "DEM", thisResDir);
    
    assembledDRGFilenameArray = writeAssembledTiles(foreDRG, backDRG, hReferenceDRG,
						    foreToReferenceOffsetDRG,  assemblerParams.weightingMode, 
						    foreToReferenceResampleFactorDRG, overwriteFlag,
						    assemblerParams.backReference,
						    pyrTileStrategyRefDRG, "DRG", thisResDir);

    //assembler - END

    
    GDALClose(hReferenceDEM);
    GDALClose(hReferenceDRG);

    for (int k = 0; k < assembledDEMFilenameArray.size(); k++){

      struct TerrainFilename assembledTerrainFilenames;
      
      if ((assembledDEMFilenameArray[k].size()!=0)&&
	  (assembledDRGFilenameArray[k].size()!=0)){
	for (int b=0; b<assembledDEMFilenameArray[k].size(); b++){
	  assembledTerrainFilenames.demFilenames.push_back(assembledDEMFilenameArray[k][b]);
	}
	for (int b=0; b<assembledDRGFilenameArray[k].size(); b++){
	  assembledTerrainFilenames.drgFilenames.push_back(assembledDRGFilenameArray[k][b]);
	}
	assembledTerrainFilenameArray.push_back(assembledTerrainFilenames);
      }
    }

    //#endif
    //#if 0
    cout<<"pyramid_assembler_antares::pyrAssembleGeoDataAntares(): numAssembledTiles="
	<<assembledTerrainFilenameArray.size()<<endl;
    
    //create the assembled tiles pyramids
    std::vector<std::vector<struct TilingParams> > pyrTileStrategyAssembledDRG;
    std::vector<std::vector<struct TilingParams> > pyrTileStrategyAssembledDEM;

    demResampleParams.pyramidMode = 0;
    drgResampleParams.pyramidMode = 0;

    for (int t = 0; t < assembledTerrainFilenameArray.size(); t++){
    
      string pyramidDirname = GetFilenameNoExt(assembledTerrainFilenameArray[t].demFilenames[0])+"_Pyr";
      cout<<"pyramidDirname="<<pyramidDirname<<endl;
            
      GDALDataset* demDS = (GDALDataset *) GDALOpen(assembledTerrainFilenameArray[t].demFilenames[0].c_str(),
						    GA_ReadOnly); 
      vector<double> tileGeographicBBox = computeGeographicBBox(demDS);
      struct bbox tileGeoBBox;
      tileGeoBBox.xl = tileGeographicBBox[0];
      tileGeoBBox.yb = tileGeographicBBox[2];
      tileGeoBBox.xr = tileGeographicBBox[1];
      tileGeoBBox.yt = tileGeographicBBox[3];
      pyrTileStrategyAssembledDEM = makeVirtualTiles(demDS, tileGeoBBox, demResampleParams);
      
      writePyramidTiles( demDS, 0,1, 
			 pyrTileStrategyAssembledDEM,
			 "DEM", pyramidDirname);
      
      int numDRGBands = 1;
      for (int b = 0; b < numDRGBands; b++){
           GDALDataset* drgDS = (GDALDataset *) GDALOpen(assembledTerrainFilenameArray[t].drgFilenames[b].c_str(),
							 GA_ReadOnly);
	   
           pyrTileStrategyAssembledDRG = changeVirtualTiles(demDS, demResampleParams, drgDS,
                                                            drgResampleParams, pyrTileStrategyAssembledDEM);
	   
	   writePyramidTiles( drgDS, 0,1, 
			      pyrTileStrategyAssembledDRG,
			      "DRG", pyramidDirname);
	   
           GDALClose((GDALDatasetH) drgDS);
      }
      
      GDALClose((GDALDatasetH) demDS);
      
    }

    //#endif
    return assembledTerrainFilenameArray;
   
  }

  //returns the set of assembled pyramids
  vector<struct TerrainFilename> pyrAssembleGeoData(string foreDEMFilename, string backDEMFilename,
						    string foreDRGFilename, string backDRGFilename,
						    struct AssemblerParams assemblerParams, 
						    string resDir, int overwriteFlag)
  {

    
    string thisResDir;
    vector<struct TerrainFilename> terrainFilenameArray;
    //vector<vector<string> >assembledTileFilenameArray;
    //assembledTileFilenameArray.resize(2);
    
    //read the fore data  
    GDALDataset* foreDEM = (GDALDataset *) GDALOpen(foreDEMFilename.c_str(), GA_ReadOnly); 
    GDALDataset* foreDRG = (GDALDataset *) GDALOpen(foreDRGFilename.c_str(), GA_ReadOnly);
    
    //compute fore bounding box in geographic coordinates
    vector<double> foreBBoxGeographic = computeGeographicBBox(foreDEM);
   
    //select overlaping tiles
    vector<string> backDEMFilenameArray = AccessDataFilesFromInput(backDEMFilename, "tif"); 
    vector<string> overlapBackDEMFilenameArray;
    overlapBackDEMFilenameArray = findOverlapFiles(backDEMFilenameArray, foreDEM);

    //compute the union bounding box in geographc coordinates
    vector<float> unionBBoxGeographic;
    vector<string> DEMFilenameArray;
    DEMFilenameArray.resize(overlapBackDEMFilenameArray.size()+1);
    DEMFilenameArray[0]=foreDEMFilename;
    for (int k=0; k < overlapBackDEMFilenameArray.size(); k++){
      DEMFilenameArray[k+1]=overlapBackDEMFilenameArray[k];
    }
    unionBBoxGeographic = computeUnionGeographicBBox(DEMFilenameArray);
    struct bbox unionGeoBBox;
    unionGeoBBox.xl = unionBBoxGeographic[0];
    unionGeoBBox.yb = unionBBoxGeographic[2];
    unionGeoBBox.xr = unionBBoxGeographic[1];
    unionGeoBBox.yt = unionBBoxGeographic[3];
    
    vector<string> overlapBackDRGFilenameArray;
    vector<string> backDRGFilenameArray = AccessDataFilesFromInput(backDRGFilename, "tif"); 
    overlapBackDRGFilenameArray = findOverlapFiles(backDRGFilenameArray, foreDRG);
    
    for (unsigned int k = 0; k < overlapBackDEMFilenameArray.size(); k++){
      
      //read the background file
      GDALDataset* backDEM = (GDALDataset *) GDALOpen(overlapBackDEMFilenameArray[k].c_str(), GA_ReadOnly); 
      GDALDataset* backDRG = (GDALDataset *) GDALOpen(overlapBackDRGFilenameArray[k].c_str(), GA_ReadOnly); 

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

      //determine the offset between the two images;
      float offsetDEM = 0;

      
      //registration
      offsetDEM = computeOffset(foreDEM, backDEM);
      string regFilename = GetFilenameNoPath(GetFilenameNoExt(foreDEMFilename))+"_reg.txt";
      saveRegistrationParams(offsetDEM, thisResDir+"/"+regFilename);
      
      
      terrainFilenameArray = pyrAssembleGeoData(foreDEM, backDEM, foreDRG, backDRG,
						unionGeoBBox, offsetDEM,
						assemblerParams, thisResDir, overwriteFlag);
      
      //save assembledFilenames to file
      string assembledDEMTileFilenameList=thisResDir+"/"+"dem_list.txt";
      ofstream foutdem(assembledDEMTileFilenameList.c_str());
      for (int i=0; i<terrainFilenameArray.size(); i++){
          foutdem << terrainFilenameArray[i].demFilenames[0] << endl;
      }
      foutdem.close();

      string assembledDRGTileFilenameList=thisResDir+"/"+"drg_list.txt";
      ofstream foutdrg(assembledDRGTileFilenameList.c_str());
      for (int i=0; i<terrainFilenameArray.size(); i++){
	for (int b=0; b<terrainFilenameArray[i].drgFilenames.size(); b++){
          foutdrg << terrainFilenameArray[i].drgFilenames[b] << endl;
	}
      }
      foutdrg.close();
  
      GDALClose(backDEM);
      GDALClose(backDRG);
    }

    GDALClose(foreDEM);
    GDALClose(foreDRG);

    return terrainFilenameArray;
  }

}





