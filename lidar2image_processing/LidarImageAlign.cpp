// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

//Eigen includes
#include <Eigen/SVD>

//GDAL includes
#include <gdal_priv.h>

//Boost include to support isisInterface
#include <boost/filesystem.hpp>

//ATK includes
#include "../common/StringUtils.h"
#include "../common/ImageProcessing.h"
#include "../camera_models/IsisInterfaceATK.h"
#include "../camera_models/PVLRead.h"
#include "../geotif_processing/CoordTransform.h"
#include "../geotif_processing/GeoUtils.h"
#include "../lidar_processing/TracksLOLA.h"
#include "../lidar_processing/TracksGCP.h"
#include "LidarImageAlign.h"

using namespace std;

namespace at
{

  void UpdateGCP( const vector<vector<AlignedLOLAShot> >& LOLATracks, 
                const string&                           camFilename, 
                vector<gcp>&                            gcpArray)
{
  int featureIndex = 0;
  int validFeatureIndex = 0;

  for( unsigned int t = 0; t < LOLATracks.size(); ++t ){
    for( unsigned int s = 0; s < LOLATracks[t].size(); ++s ){

      if( LOLATracks[t][s].featurePtLOLA == 1 ){

	gcpArray[featureIndex].filename.push_back(camFilename);
	
        if( (LOLATracks[t][s].valid == 1) && 
            (LOLATracks[t][s].reflectance != 0) && 
            (LOLATracks[t][s].reflectance != -1)  ){

          // convert a map projected pixel location to the
          // original image coordinate system.
          gcpArray[featureIndex].x.push_back( LOLATracks[t][s].image_x);
          gcpArray[featureIndex].y.push_back( LOLATracks[t][s].image_y);

          gcpArray[featureIndex].x_before.push_back(LOLATracks[t][s].imgPt[2].x);
          gcpArray[featureIndex].y_before.push_back(LOLATracks[t][s].imgPt[2].y);

          gcpArray[featureIndex].trackIndex = t;
          gcpArray[featureIndex].shotIndex = s;

          validFeatureIndex++;
        }//valid==1
	else
	{
	  gcpArray[featureIndex].x.push_back(-1);
	  gcpArray[featureIndex].y.push_back(-1);
	  gcpArray[featureIndex].x_before.push_back(-1);
	  gcpArray[featureIndex].y_before.push_back(-1);
	  //gcpArray[featureIndex].trackIndex = -1;
	  //gcpArray[featureIndex].shotIndex = -1;
	}
        featureIndex++;
      }
    }
  }
}


  //determines all points in all tracks and all shots that overlap with
  //the cub image
  int GetAllPtsFromCub( vector<vector<LOLAShot> >& LOLATracks,
		      IsisInterface* cam, 
		      GDALDataset* camData, struct bbox imageROI)
{
    int numValidImgPts = 0;

    vector<double> xyz; xyz.resize(3);
    
    int width  = GDALGetRasterXSize(camData);
    int height = GDALGetRasterYSize(camData);
    GDALRasterBand* camBandData = camData->GetRasterBand(1);
    float *image = new float[height*width];
    camBandData->RasterIO(GF_Read, 0, 0, width, height,  
			  image, width, height, 
			  GDT_Float32, 0, 0);
    
    for( unsigned int k = 0; k < LOLATracks.size(); ++k ){
      for( unsigned int i = 0; i < LOLATracks[k].size(); ++i ){
	vector<LOLAPoint> points = LOLATracks[k][i].LOLAPts;
	
	LOLATracks[k][i].image = 1; 
	LOLATracks[k][i].imgPt.resize( points.size() );
	
	for( unsigned int j = 0; j < points.size(); ++j ){
	  
	  Eigen::Vector3f lonLatRad;
	  lonLatRad[0] = M_PI*points[j].x/180.0;
	  lonLatRad[1] = M_PI*points[j].y/180.0;
	  lonLatRad[2] = 1000.0*points[j].z;

	  Eigen::Vector3f xyz_temp = geographicToCartesian(lonLatRad);
	  
	  xyz[0]=xyz_temp[0];
	  xyz[1]=xyz_temp[1];
	  xyz[2]=xyz_temp[2];
	  
	  vector<double> cub_pix = cam->point_to_pixel(xyz);
	  if ((cub_pix[0]==0)||(cub_pix[1]==0)){
	    cout<<"invalid pixel!"<<endl;
	  }		 
	  if ((cub_pix[0]>imageROI.xl) && (cub_pix[0]<imageROI.xr) &&
	      (cub_pix[1]>imageROI.yt) && (cub_pix[1]<imageROI.yb)){ 
            LOLATracks[k][i].imgPt[j].val = GetInterpValue(cub_pix[0], cub_pix[1], image, width, height);
	    LOLATracks[k][i].imgPt[j].x = cub_pix[0];
	    LOLATracks[k][i].imgPt[j].y = cub_pix[1];

	    numValidImgPts++;
	  }
	  else{//invalidate the point
	    LOLATracks[k][i].image = -1;
	    break;
	  }
	}
      }//i  
    }//k

    delete image;
    
    cout<<"numValidImgPts="<<numValidImgPts<<endl;
    return numValidImgPts;
  }

  //make a complete copy of all tracks (valid and invalid)
  //fill in the AlignedLOLAShot specific fields (image_x, image_y and image)
  std::vector<std::vector< AlignedLOLAShot> > initialize_aligned_lola_shots(std::vector<std::vector<LOLAShot> >& tracks)
  {
    std::vector<std::vector< AlignedLOLAShot> > alignedTracks;

    for (unsigned int i = 0; i < tracks.size(); i++){
      std::vector< AlignedLOLAShot> alignedTrack;
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	AlignedLOLAShot alignedShot(tracks[i][j]);
        if (!tracks[i][j].valid || tracks[i][j].imgPt.size() <= 2){
	    alignedShot.image_x = -1;
	    alignedShot.image_y = -1;
	    alignedShot.image   = -1;
            //alignedShot.valid   = -1;
	}
	else{
	    alignedShot.image_x = tracks[i][j].imgPt[2].x;
	    alignedShot.image_y = tracks[i][j].imgPt[2].y;
	    alignedShot.image   = tracks[i][j].imgPt[2].val;
            //alignedShot.valid = 1;
	}
	alignedTrack.push_back(alignedShot);
      }
      alignedTracks.push_back(alignedTrack);
    }

    return alignedTracks;
  }

  int computeNumberOfValidPoints(vector<AlignedLOLAShot> & track)
  {

    int numValidPoints = 0;
    for (unsigned int i = 0; i < track.size(); i++){
      if ((track[i].image != -1) && (track[i].valid == 1)){	
        numValidPoints++;    
      }
    }

    return numValidPoints;
  }

  Eigen::MatrixXd computeReconstructionErrorVec(vector<AlignedLOLAShot> & track, int numValidPoints)
  {

    Eigen::MatrixXd errorVec(numValidPoints, 1);
    int validPtIndex = 0;

    for (unsigned int i = 0; i < track.size(); i++){

      if ((track[i].image != -1) && (track[i].valid == 1)){	
	errorVec(validPtIndex, 0) = track[i].image - track[i].synth_image;

        if (isnan(errorVec(validPtIndex))){
          cout<<"LidarImageAlign.cc::computeReconstructionErrorVec(): ERROR!"<<endl;
          exit(1);
	}
	if (fabs(errorVec(validPtIndex, 0)) > 1.0){
          cout<<"Warning!!!! errorVec["<<validPtIndex<<"]="<<errorVec(validPtIndex, 0) <<endl;
	  //exit(1);
	  errorVec(validPtIndex, 0) = 1.0;
	}
        validPtIndex++;
      }
       
    }
    if (validPtIndex != numValidPoints){
      cout<<"LidarImageAlign.cc::computeReconstructionErrorVec(): ERROR! validPtIndex="<<validPtIndex
          <<", numValidPoints="<<numValidPoints<<endl;
      exit(1);
    }
    return errorVec;
  }

  float computeReconstructionErrorLS(Eigen::MatrixXd errorVec)
  {
    float errorLS = 0.0;
    for (unsigned int i = 0; i < errorVec.rows(); i++){
      if (isnan(errorVec(i,0))){
	  cout<<"LidarImageAlign.cc::computeReconstructionErrorLS(): Error! in element "<<i<<endl;
	  exit(1);
      }
      errorLS = errorLS + errorVec(i, 0)*errorVec(i, 0);  
    }
    if (isnan(errorLS)){
      cout<<"LidarImageAlign.cc::computeReconstructionErrorLS(): Error!"<<endl;
      exit(1);
    }
    return errorLS;
  }

  Eigen::MatrixXd computeJacobianAffine(vector<AlignedLOLAShot> & track,
				       float* image, int width, int height, 
                                       int numPoints)
  {
    Eigen::MatrixXd J(numPoints, 6); 
    double h = 1.0;//this should be in the configuration file

    int index = 0;

    for (unsigned int i = 0; i < track.size(); i++){

      if ((track[i].image !=-1) && (track[i].valid != 0)){
      
	int x = track[i].image_x;
	int y = track[i].image_y;
		
	//fininte difference with three points to right and left
	double dx, dy;
	if (x >= width - 3 || x <= 3 || y >= height - 3 || y <= 3){
	  dx = 0.0; dy = 0.0;
          cout<<"computeJacobianAffine(): I should not be here!"<<endl;
          exit(1);
	}
	else{
	  
	  dx = (-1.0/60)* image[y*width+(x-3)]
	     + (3.0/20) * image[y*width+(x-2)]
	     - (3.0/4)  * image[y*width+(x-1)]
	     + (3.0/4)  * image[y*width+(x+1)]
	     - (3.0/20) * image[y*width+(x+2)]
	     + (1.0/60) * image[y*width+(x+3)];
	  
	  dy = (-1.0/60)* image[(y-3)*width + x]
	     + (3.0/20) * image[(y-2)*width + x]
	     - (3.0/4)  * image[(y-1)*width + x]
	     + (3.0/4)  * image[(y+1)*width + x]
	     - (3.0/20) * image[(y+2)*width + x]
	     + (1.0/60) * image[(y+3)*width + x];
	  
	}
	
	//compute the derivatives
	J(index, 0) = dx / (h / track[i].imgPt[2].x);
	J(index, 1) = dx / (h / track[i].imgPt[2].y);
	J(index, 2) = dx /  h;
	J(index, 3) = dy / (h / track[i].imgPt[2].x);
	J(index, 4) = dy / (h / track[i].imgPt[2].y);
	J(index, 5) = dy /  h;
	index++;
      }
    }
    return J;
  }


  //Gauss Newton solution
  //returns the matching error r_err; 
  Eigen::Matrix3d gaussNewtonAffine(vector<AlignedLOLAShot> & track,
				    float* image, int width, int height,
				    struct bbox imageROI, Eigen::Matrix3d matrix,                                      
				    float *matchingError)
  {
    Eigen::Matrix3d B = matrix;
    cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): B="<<endl;
    cout<<B<<endl;

    int numValidPoints = computeNumberOfValidPoints(track);
    Eigen::MatrixXd errorVec = computeReconstructionErrorVec(track, numValidPoints);
    float currError = computeReconstructionErrorLS(errorVec);
    float prevError = currError;

    cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): iteration=0, numInitValidPoints="<<numValidPoints
        <<", currError="<<currError<<endl;
    
    if (numValidPoints <= 0){ // no points
      *matchingError = currError;
      return B;
    }
     
 
    for (int iteration = 0; iteration < MAX_GAUSS_NEWTON_STEPS; iteration++){
     
    
      Eigen::MatrixXd J  = computeJacobianAffine(track, image, width,
						 height, numValidPoints);
      Eigen::MatrixXd JT = J.transpose();
      
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(JT*J, Eigen::ComputeThinU | Eigen::ComputeThinV);
      const Eigen::VectorXd s  = svd.singularValues();
      const Eigen::MatrixXd U  = svd.matrixU();
      const Eigen::MatrixXd V  = svd.matrixV();
      
      Eigen::MatrixXd S(6, 6);
      for (int i = 0; i < 6; i++){
	for (int j = 0; j < 6; j++){
	  S(i,j)=0;
	}
      }
      
      for (int i = 0; i < 6; i++){
	if (s(i) > 10e-15){ // don't divide by 0
	  S(i, i) = 1.0 / s(i);
	}
      }

      //cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): S="<<endl;
      //cout<<S<<endl;

      Eigen::MatrixXd pseudoinverse = V * S * U.transpose();
      
      Eigen::MatrixXd delta = - pseudoinverse * (JT * errorVec);

      Eigen::Matrix3d prevMatrix = B;

      if ( isnan(delta(0,0)) || isnan(delta(1,0)) || isnan(delta(2,0)) || isnan(delta(3,0)) || isnan(delta(4,0)) || isnan(delta(5,0)) ){
          cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): delta="<<endl;
	   *matchingError = prevError;  
	   return prevMatrix;

      }    
      //cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): delta="<<endl;
      //cout<<delta<<endl;

     

      B(0, 0) += delta(0, 0); B(0, 1) += delta(1, 0); B(0, 2) += delta(2, 0);
      B(1, 0) += delta(3, 0); B(1, 1) += delta(4, 0); B(1, 2) += delta(5, 0);

      //cout<<"LidarImageAlign.cc:: gaussNewtonAffine(): B="<<endl;
      //cout<<B<<endl;
      
  
      update_aligned_track(track, B, image, width, height, imageROI);
      
      numValidPoints = computeNumberOfValidPoints(track);
      errorVec = computeReconstructionErrorVec(track, numValidPoints);
      currError = computeReconstructionErrorLS(errorVec);
     
      cout<<"LidarImageAlign.cc:: gaussNewtonAffine: iteration = "<< iteration<<", error = "<<currError<<endl;
      //cout<<"LidarImageAlign.cc:: gaussNewtonAffine: err="<<prevError<<" before iter "<<iteration<<endl;
      //cout<<"LidarImageAlign.cc:: gaussNewtonAffine: err ="<<currError<<" after iter  "<<iteration<<endl;
      
      if (currError == -1){
	cout<<"LidarImageAlign.cc:: gaussNewtonAffine: Matrix jumped too far, exiting after "<<iteration<<" iterations."<<endl;
        *matchingError = prevError;  
	return prevMatrix;
      }
      
      if (currError > 1.004*prevError){
	cout<<"LidarImageAlign.cc:: gaussNewtonAffine: Errors are increasing, exiting after "<<iteration<<" iterations."<<endl;
	*matchingError = prevError;  
	return prevMatrix;
      }
      
      /*
      if ((fabs((currError-prevError)/numValidPoints) < 0.000001) && (iteration>1)){
	*matchingError = currError;  
	return B;
      }
      */

      prevError = currError;

    }

    *matchingError = currError;
   
     return B;

  }

  vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> >& tracks,
						  GDALDataset *imageDS, struct bbox imageROI,
						  Eigen::Matrix3d& trans, float *error)
  {
 
    int imageWidth  = GDALGetRasterXSize(imageDS);
    int imageHeight = GDALGetRasterYSize(imageDS);
    GDALRasterBand* imageBand = imageDS->GetRasterBand(1);
    float *imageData = new float[imageWidth*imageHeight];
    
    imageBand->RasterIO(GF_Read, 0, 0, imageWidth, imageHeight,  
			imageData, imageWidth, imageHeight, 
			GDT_Float32, 0, 0);
    
    vector<vector< AlignedLOLAShot> > alignedTracks = initialize_aligned_lola_shots(tracks);
    
    cout<<"LidarImageAlign.cc:: align_to_image(): transform before align="<<endl;
    cout<<trans<<endl;
    
    update_aligned_tracks(alignedTracks, trans, imageData, imageWidth, imageHeight, imageROI);
    //cout<<"LidarImageAlign.cc:: align_to_image(): done transform_tracks"<<endl;

    //float initError = compute_reconstruct_error_ls(alignedTracks);
    
    std::cout.precision(12);
   
    cout<<"LidarImageAlign.cc:: align_to_image(): Append tracks"<<endl;
    vector<AlignedLOLAShot> oneAlignedTrack;
    for (unsigned int i = 0; i < alignedTracks.size(); i++){
      oneAlignedTrack.insert(oneAlignedTrack.end(), alignedTracks[i].begin(), alignedTracks[i].end());
    }

    int numValidPoints = computeNumberOfValidPoints(oneAlignedTrack);
    Eigen::MatrixXd errorVec = computeReconstructionErrorVec(oneAlignedTrack, numValidPoints);
    float initError = computeReconstructionErrorLS(errorVec);
  
    cout<<"-----------------------------------------------------"<<endl;
    cout<<"LidarImageAlign.cc:: align_to_image(): Global Gauss-Newton alignment"<<endl;
    cout<<"LidarImageAlign.cc:: align_to_image(): transform before align="<<endl;
    cout<<trans<<endl;
    cout<<"LidarImageAlign.cc:: align_to_image(): initial error = "<<initError<<endl;

    //trans = gauss_newton_affine(oneAlignedTrack, imageData, imageWidth, imageHeight,
    //				imageROI, trans, error);
    float finalError;
    trans = gaussNewtonAffine(oneAlignedTrack, imageData, imageWidth, imageHeight,
    			      imageROI, trans, &finalError);

    cout<<"LidarImageAlign.cc:: align_to_image(): transform after align="<<endl;
    cout<<trans<<endl;
    cout<<"LidarImageAlign.cc:: align_to_image(): final error = "<<finalError<<endl;
    cout<<"-----------------------------------------------------"<<endl;
    if (finalError < initError){
      *error = finalError;
    }
    else{
      *error = initError;
    }
    
    //updates the set of aligned lola tracks (multiple tracks)
    //updates the reconstructed image for a linear factor (gain)
    //updates shot.image_x, shot.image_y and shot.image
    //shot.image=-1 for invalid points (missing reflectance or outside image boundaries)
    update_aligned_tracks(alignedTracks, trans, imageData, imageWidth, imageHeight, imageROI);
    
    delete imageData;
    
    return alignedTracks;
  }


  
  vector<vector<AlignedLOLAShot> > align_to_image_pyramid(vector<vector<LOLAShot> > & trackPts, 
                                                          const string & imageFilename,
                                                          Eigen::Matrix3d & trans,
                                                          int *r_numValidReflectancePoints,
                                                          float *error,
							  string resultDirname)
  {

    int zoomMultiplier = 2;  
    int zoomFactor = 16; //this should go to the settings file.
    
    int initZoomFactor = zoomFactor;
    
    vector<vector< AlignedLOLAShot> > aligned;
    
    bool first = true;
    
    GDALAllRegister();
    
    GDALDataset* camData = (GDALDataset *) GDALOpen(imageFilename.c_str(), GA_ReadOnly);
    
    //read the cub image.
    int camWidth  = GDALGetRasterXSize(camData);
    int camHeight = GDALGetRasterYSize(camData);
    GDALRasterBand* camBandData = camData->GetRasterBand(1);
    GDALDataType typeCamData  = camBandData->GetRasterDataType();
    cout<<"LidarImageAlign.cc:: cam width="<<camWidth<<", height="<<camHeight<<endl;

    short int *camImageC16 = NULL;
    float     *camImageF32 = NULL;
    float      minImageVal;
    float      maxImageVal;

    GDALDataset *camData32F;
    
    if (typeCamData==GDT_Int16){
      string image32FFilename;
      if (resultDirname.compare("")==0){
	image32FFilename = GetFilenameNoExt(GetFilenameNoPath(imageFilename))+"_F32.tiff";
	camData32F = uint16ToF32GDAL(camData, image32FFilename, 0);
      }
      else{
	image32FFilename = resultDirname+"/"+GetFilenameNoExt(GetFilenameNoPath(imageFilename))+"_F32.tiff";
	camData32F = uint16ToF32GDAL(camData, image32FFilename, 1);
      }    
      GDALClose((GDALDatasetH) camData);
    }
    if (typeCamData==GDT_Float32){
      GDALClose((GDALDatasetH) camData);
      camData32F = (GDALDataset *) GDALOpen(imageFilename.c_str(), GA_ReadOnly);
    }
    
    GDALRasterBand* camBandData32F = camData32F->GetRasterBand(1);
 
    float noDataValue = camBandData32F->GetNoDataValue();
   
    //compute min max image range
    minImageVal = camBandData32F->GetMinimum();
    maxImageVal = camBandData32F->GetMaximum();
     
    cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): minImageVal="<<minImageVal<<", maxImageVal="<<maxImageVal<<endl;

    string pvlFilename=GetFilenameNoExt(imageFilename)+".pvl";
    cout<<"pvlFilename="<<pvlFilename<<endl;
    Eigen::Vector3f lightPosition =scanPVLForSunPosition(pvlFilename);
    Eigen::Vector3f cameraPosition=scanPVLForCamPosition(pvlFilename);

    cout<<"before opening "<<imageFilename<<endl;
    IsisInterface *model = IsisInterface::open(imageFilename); 
    cout<<"after opening "<<imageFilename<<endl;
    struct bbox imageROI;
    imageROI.xl = 3*zoomFactor;
    imageROI.yt = 3*zoomFactor;
    imageROI.xr = imageROI.xl + camWidth-2*3*zoomFactor;
    imageROI.yb = imageROI.yt + camHeight-2*3*zoomFactor;
    
    int numValidImgPoints = GetAllPtsFromCub(trackPts, model, camData32F, imageROI);
    delete model;
     
    int numValidReflectancePoints = ComputeAllReflectance(trackPts, cameraPosition, lightPosition);
    *r_numValidReflectancePoints = numValidReflectancePoints;

    GDALDataset* displayData32F;
    GDALDriver*  displayDriver;

    if (resultDirname.compare("")!=0){
      
      //draw the initial tracks (only those within the image boundaries)
      displayDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
      string initTracksImageFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))+"_init_tracks.jpg";
      
      displayData32F = displayDriver->CreateCopy(initTracksImageFilename.c_str(), camData32F,
						 FALSE, NULL, NULL, NULL);
      DrawTracksOnImage(trackPts, displayData32F, "tracks", 8);
      GDALClose(displayData32F);
            
      //draw the initial reflectance (only for points with valid reflectance)
      string InitReflTracksImageFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))+
	"_init_refl_tracks.jpg";
      displayData32F = displayDriver->CreateCopy(InitReflTracksImageFilename.c_str(), camData32F,
						 FALSE, NULL, NULL, NULL);
      DrawTracksOnImage(trackPts, displayData32F, "reflectance", 8);
      GDALClose(displayData32F);
      //draw the initial tracks - END
    }
  
    Eigen::Matrix3d thisMat;
    thisMat(0,0) = 1.0/zoomFactor;  thisMat(0,1) = 0;                  thisMat(0,2) = 0;
    thisMat(1,0) = 0;               thisMat(1,1) = 1.0/zoomFactor;     thisMat(1,2) = 0;
    thisMat(2,0) = 0;               thisMat(2,1) = 0;                  thisMat(2,2) = 0;
    
    cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): zoomFactor="<<zoomFactor<<endl;
    transform_tracks_by_matrix(trackPts, thisMat, imageROI);
    
    while (zoomFactor >= 1.0){
      
      if (first){
	first = false;
      }
      else{
        Eigen::Matrix3d t1;
	t1(0,0) = (float)zoomMultiplier; t1(0,1) = 0;                        t1(0,2) = 0;
	t1(1,0) = 0;                      t1(1,1) = (float)zoomMultiplier;   t1(1,2) = 0;
	t1(2,0) = 0;                      t1(2,1) = 0;                       t1(2,2) = 1;
	Eigen::Matrix3d t2;
	t2(0,0) = 1.0/zoomMultiplier;     t2(0,1) = 0;                       t2(0,2) = 0;
	t2(1,0) = 0;                      t2(1,1) = 1.0/zoomMultiplier;      t2(1,2) = 0;
	t2(2,0) = 0;                      t2(2,1) = 0;                       t2(2,2) = 1;

	cout<<"zoomFactor="<<zoomFactor<<", new trans="<<endl;
	cout<<trans<<endl;
	
	trans = t1 * trans * t2;
	cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): zoomFactor="<<zoomFactor<<", new trans="<<endl;
	cout<<trans<<endl;
	
	Eigen::Matrix3d zoomMatrix;
	zoomMatrix(0,0) = zoomMultiplier;  zoomMatrix(0,1) = 0;               zoomMatrix(0,2) = 0;
	zoomMatrix(1,0) = 0;               zoomMatrix(1,1) = zoomMultiplier;  zoomMatrix(1,2) = 0;
	zoomMatrix(2,0) = 0;               zoomMatrix(2,1) = 0;               zoomMatrix(2,2) = 1;
	
        cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): zoomFactor="<<zoomFactor<<", tempM="<<endl;
	cout<<zoomMatrix<<endl;
	
	transform_tracks_by_matrix(trackPts, zoomMatrix, imageROI);
      }
      
      cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): zoomFactor="<<zoomFactor<<endl;
      
      GDALDataset* zoomImageDS;
      string filename;
      if (resultDirname.compare("")!=0){
          filename = resultDirname+"/temp_"+NumberToString((int)floor(zoomFactor))+".tif";
      }
      else{
	  filename = "temp_"+NumberToString((int)floor(zoomFactor))+".tif";
      }
      GDALDatasetH hZoomImage = resampleGDALDataset(camData32F, 1.0/zoomFactor, filename, string("floor"), 0);
     
      zoomImageDS=(GDALDataset*)hZoomImage;
      GDALRasterBandH hZoomImageBand = GDALGetRasterBand(zoomImageDS, 1);
      GDALRasterBand *zoomImageBand  = zoomImageDS->GetRasterBand(1);
      zoomImageBand->SetNoDataValue(camBandData32F->GetNoDataValue());
      int zoomImageWidth  = GDALGetRasterXSize(zoomImageDS);
      int zoomImageHeight = GDALGetRasterYSize(zoomImageDS);
      CPLErr imError   = CE_None;
      float*zoomImage32F = new float[zoomImageWidth * zoomImageHeight];
      
      camBandData32F->RasterIO(GF_Read, 0, 0, camWidth, camHeight, zoomImage32F,
			       zoomImageWidth, zoomImageHeight, GDT_Float32, 0, 0);
      imError = GDALRasterIO(hZoomImageBand, GF_Write, 0, 0, zoomImageWidth, zoomImageHeight,
			     zoomImage32F, zoomImageWidth, zoomImageHeight , GDT_Float32, 0, 0);
      
      cout<<"gdal: "<<zoomImageWidth<<", "<<zoomImageHeight<<endl;
      cout<<"cv:   "<<camWidth/zoomFactor<<", "<<camHeight/zoomFactor<<endl;
        
      //struct bbox imageROI;
      imageROI.xl = 3*zoomFactor; 
      imageROI.yt = 3*zoomFactor;
      imageROI.xr = imageROI.xl + zoomImageWidth-2*3*zoomFactor;
      imageROI.yb = imageROI.yt + zoomImageHeight-2*3*zoomFactor;
      
      //original trackPts zoomed only, current affine transformation matrix
      aligned = align_to_image(trackPts, zoomImageDS, imageROI, trans, error);
  
      //drawing - START
      
      if (resultDirname.compare("")!=0){

	GDALDataset* displayZoomImageDS;
	
	//draw reconstructed image at each level of the pyramid
	string synthImageTracksFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))
	  + "_aligned_refl_"+NumberToString((int)floor(zoomFactor))+".jpg";
	displayZoomImageDS = displayDriver->CreateCopy(synthImageTracksFilename.c_str(), zoomImageDS,
						       FALSE, NULL, NULL, NULL);
	
	int thickness = 8/zoomFactor; if (thickness==0){thickness=1;}
	DrawAlignedImageTracks(aligned, displayZoomImageDS, "reflectance", thickness);
	GDALClose(displayZoomImageDS);
	
	//draw tracks-only at each level of the pyramid
	string imageTracksFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))
	  + "_aligned_tracks_"+NumberToString((int)floor(zoomFactor))+".jpg";
	displayZoomImageDS = displayDriver->CreateCopy(imageTracksFilename.c_str(), zoomImageDS,
						       FALSE, NULL, NULL, NULL);
	DrawAlignedImageTracks(aligned, displayZoomImageDS, "tracks", thickness);
	GDALClose(displayZoomImageDS);
	//drawing - START
      }

      cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): zoom_factor="<<zoomFactor<<endl;
      cout<<"LidarImageAlign.cc:: align_to_image_pyramid(): outTrans="<<endl;
      cout<<trans<<endl;
	
      zoomFactor = zoomFactor/zoomMultiplier;
      
      GDALClose(zoomImageDS);

    }
  
    return aligned;
  }




  void printRHS_LHS( Eigen::MatrixXf rhs,  Eigen::VectorXf lhs, int index, int iter)
  {
    cout<<"index= "<<index<<", iter= "<<iter<<endl;
    printf("--------------------------------------\n");
    printf("RHS = \n");
    printf("[ %f %f %f | %f %f %f]\n",rhs(0,0), rhs(0,1), rhs(0,2), rhs(0,3), rhs(0,4), rhs(0,5) );   
    printf("[ %f %f %f | %f %f %f]\n", rhs(1,0), rhs(1,1), rhs(1,2), rhs(1,3), rhs(1,4),rhs(1,5) ); 
    printf("[ %f %f %f | %f %f %f]\n", rhs(2,0), rhs(2,1),rhs(2,2), rhs(2,3), rhs(2,4), rhs(2,5) );
    printf("[ %f %f %f | %f %f %f]\n", rhs(3,0), rhs(3,1), rhs(3,2), rhs(3,3), rhs(3,4), rhs(3,5) );
    printf("[ %f %f %f | %f %f %f]\n", rhs(4,0), rhs(4,1), rhs(4,2), rhs(4,3), rhs(4,4), rhs(4,5) );
    printf("[ %f %f %f | %f %f %f]\n", rhs(5,0), rhs(5,1), rhs(5,2), rhs(5,3), rhs(5,4), rhs(5,5) );
    printf("--------------------------------------\n\n");
    printf("LHS = [ %f, %f, %f, %f, %f, %f]\n", lhs[0], lhs[1], lhs[2], lhs[3], lhs[4], lhs[5]);
    printf("--------------------------------------\n\n");
  }

  void printIterationValues(  float error, int index, int iter, int numValidPts)
  {
    cout<<"index= "<<index<<", iter= "<<iter<<", error="<<error<<", numValidPts="<<numValidPts<<endl;
  }
  /*
    //This functions needs to be changed to generate multiple starting points for each 
    //element of the affine transform 
    void GenerateInitTransforms( vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
    {
    initTransfArray.resize(settings.maxNumStarts);
    for (int i = 0; i < settings.maxNumStarts; i++){
    initTransfArray[i][0] = 1.0;
    initTransfArray[i][1] = 0.0;
    initTransfArray[i][2] = (i-settings.maxNumStarts/2)*2;
    initTransfArray[i][3] = 0.0;
    initTransfArray[i][4] = 1.0;
    initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }  
    }
  */

  void GainBiasAccumulator( const vector<AlignedLOLAShot>& trackPts, 
			    float&            sum_rfl, 
			    float&            sum_img, 
			    float&            sum_rfl_2, 
			    float&            sum_rfl_img, 
			    int&              numValidPts )
  {
    
    for (unsigned int i = 0; i < trackPts.size(); ++i){    
      if ((trackPts[i].reflectance >= 0) && trackPts[i].image >= 0){
	//update the nominator for the center point
	sum_rfl     += trackPts[i].reflectance;
	sum_rfl_2   += trackPts[i].reflectance*trackPts[i].reflectance;
	sum_rfl_img += trackPts[i].reflectance*(trackPts[i].image);
	sum_img     += trackPts[i].image;
	++numValidPts;
      }
    }
    
  }
  
  Eigen::Vector2f GainBiasSolver( const float& sum_rfl, 
                        const float& sum_img, 
                        const float& sum_rfl_2, 
                        const float& sum_rfl_img, 
                        const int&   numValidPts ){
  /*Matrix<float,2,2> rhs;
  Vector<float,2> lhs;

  rhs(0,0) = sum_rfl_2;
  rhs(0,1) = sum_rfl;
  rhs(1,0) = sum_rfl;
  rhs(1,1) = numValidPts;
  lhs(0) = sum_rfl_img;
  lhs(1) = sum_img;
  solve_symmetric_nocopy(rhs,lhs);
  return lhs;*/
  //Vector<float,2> 
  Eigen::Vector2f result;
  result(0) = sum_img / sum_rfl;
  result(1) = 0;
  return result;
}


  //computes the gain and bias factor for each track
  //Vector2 
  Eigen::Vector2f ComputeGainBiasFactor( const vector<AlignedLOLAShot>& trackPts )
  {
    int numValidPts = 0;
    float sum_rfl = 0.0; 
    float sum_img = 0.0;
    float sum_rfl_2 = 0.0;
    float sum_rfl_img = 0.0;
    //Vector2 gain_bias;
    Eigen::Vector2f gain_bias;
    GainBiasAccumulator( trackPts, sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
    
    //cout<<"NUM_VALID_POINTS="<<numValidPts<<endl;
    //if (numValidPts != 0){
    if( numValidPts > 1 ){ 
      gain_bias = GainBiasSolver( sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
    }
    else{
      //invalid scaleFactor, all tracks are invalid
      gain_bias(0) = 0.0;
      gain_bias(1) = 0.0;
    }
    return gain_bias;
  }

  //computes the gain and bias factor for all tracks at once
  Eigen::Vector2f ComputeGainBiasFactor( const vector<vector<AlignedLOLAShot> >& trackPts )
  {
    
    cout<<"ComputeGainBiasFactor"<<endl;  
    int numValidPts = 0;
    float sum_rfl = 0.0; 
    float sum_img = 0.0;
    float sum_rfl_2 = 0.0;
    float sum_rfl_img = 0.0;
    Eigen::Vector2f gain_bias;
    
    for (unsigned int k = 0; k < trackPts.size(); ++k){
      GainBiasAccumulator( trackPts[k], sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
    }
    
    if (numValidPts != 0){ 
      gain_bias = GainBiasSolver( sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
      
    }
    else{
      //invalid scaleFactor, all tracks are invalid
      gain_bias(0) = 0.0;
      gain_bias(1) = 0.0;
    }
    
    return gain_bias;
  }
  
  void apply_gain_bias(vector<vector<AlignedLOLAShot> >& tracks)
  {
    Eigen::Vector2f gain_bias = ComputeGainBiasFactor( tracks );
    cout<<"LidarImageAlign.cc:: apply_gain_bias(multi-track): gain="<<gain_bias(0)<<", bias="<<gain_bias(1)<<endl;

    int numPoints = 0;
    int numInvalidPoints = 0;
    for (unsigned int i = 0; i < tracks.size(); i++){
      for (unsigned int j = 0; j < tracks[i].size(); j++){

	if (tracks[i][j].valid==1){
	  
	  tracks[i][j].synth_image = tracks[i][j].reflectance * gain_bias(0) + gain_bias(1);
	  
	  numPoints++;
	  
	  if (tracks[i][j].synth_image > 1.0){
	    //cout<<"LidarImageAlign.cc:: apply_gain_bias(multi-track): LARGE IMAGE RECONSTR VAL = "
	    //	<<tracks[i][j].synth_image<<endl;
	    tracks[i][j].synth_image = 1.0;
	    numInvalidPoints++;
	  }
	  
	  if (tracks[i][j].synth_image < 0){
	    cout<<"LidarImageAlign.cc:: apply_gain_bias(multi-track): NEGATIVE IMAGE RECONSTR VAL = "
		<<tracks[i][j].synth_image<<endl;
	    tracks[i][j].synth_image = 0;
	  }
	}
      }
    }
    cout<<"apply_gain_bias DONE"<<endl;
    cout<<"numPoints="<<numPoints<<", numInvalidPoints="<<numInvalidPoints<<endl;
  }
  
  void apply_gain_bias(vector<AlignedLOLAShot>& track)
  {
    Eigen::Vector2f gain_bias = ComputeGainBiasFactor( track );
    //cout<<"LidarImageAlign.cc:: apply_gain_bias(single-track): gain="
    //	<<gain_bias(0)<<", bias="<<gain_bias(1)<<endl;

    for (unsigned int i = 0; i < track.size(); i++){

      if (track[i].valid == 1){
	track[i].synth_image = track[i].reflectance * gain_bias(0) + gain_bias(1);
	
	if (track[i].synth_image > 1.0){ 
	  track[i].synth_image = 1.0;
	}
	
	if (track[i].synth_image < 0){
	  cout<<"LidarImageAlign.cc:: apply_gain_bias(single-track): NEGATIVE IMAGE RECONSTR VAL = "
	      <<track[i].synth_image<<endl;
	  track[i].synth_image = 0;	
	}
      }
    } 
  }

  //shifts the *original* image points of all tracks
  void transform_tracks_by_matrices(vector<vector<LOLAShot> > & LOLATracks,
				    vector<Eigen::Matrix3d> matrices)
  {
    for (unsigned int i = 0; i < LOLATracks.size(); i++){
      //Matrix3x3 M = matrices[i];
      Eigen::Matrix3d M = matrices[i];
      for (unsigned int j = 0; j < LOLATracks[i].size(); j++){
	for (unsigned int k = 0; k < LOLATracks[i][j].imgPt.size(); k++){
	  Eigen::Vector3d r = M * Eigen::Vector3d(LOLATracks[i][j].imgPt[k].x, LOLATracks[i][j].imgPt[k].y, 1);
	  LOLATracks[i][j].imgPt[k].x = (float)r(0);
	  LOLATracks[i][j].imgPt[k].y = (float)r(1);
	}
      }
    }
  }
  
  //shifts the *original* image points of all the tracks
  void transform_tracks_by_matrix(vector<vector<LOLAShot> > & tracks,
				  Eigen::Matrix3d transform,
                                  bbox imageROI)
  {
    for (unsigned int i = 0; i < tracks.size(); i++){
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	for (unsigned int k = 0; k < tracks[i][j].imgPt.size(); k++){
	  Eigen::Vector3d r = transform * Eigen::Vector3d(tracks[i][j].imgPt[k].x, tracks[i][j].imgPt[k].y, 1);
	  tracks[i][j].imgPt[k].x = (float)(r(0));
	  tracks[i][j].imgPt[k].y = (float)(r(1));
          if ((tracks[i][j].imgPt[k].x < imageROI.xl)||(tracks[i][j].imgPt[k].y < imageROI.yt) ||
              (tracks[i][j].imgPt[k].x > imageROI.xr)||(tracks[i][j].imgPt[k].y > imageROI.yb)){
	    tracks[i][j].image = -1;
	    //cout<<"LidarImageAlign.cc:: transform_tracks_by_matrix(): ERROR, x="
	    //<<tracks[i][j].imgPt[k].x<<", y="<<tracks[i][j].imgPt[k].y<<endl;
	  }
	}
      }
    }
    //cout<<"exiting transform_tracks_by_matrix()"<<endl;
    //exit(1);
  }


  //adjust the pixel position associated with each point in the aligned Lidar track
  void transform_aligned_track_by_matrix(vector<AlignedLOLAShot> & track,
					 Eigen::Matrix3d transform,
					 float* image, int width, int height,
					 bbox imageROI)
  {

    for (unsigned int i = 0; i < track.size(); i++){

      if (!track[i].valid){
	track[i].image = -1;       
      }
      else{

	//apply transformation on the original points
	Eigen::Vector3d r = transform * Eigen::Vector3d(track[i].imgPt[2].x, track[i].imgPt[2].y, 1);
	r = r / r(2);
	int x = (int)r(0);
	int y = (int)r(1);
	track[i].image_x = x;
	track[i].image_y = y;

	//if (x < roi.x || y < roi.y || x >= roi.x+roi.width || y >= roi.y+roi.height){
	if (x < imageROI.xl || y < imageROI.yt || x >= imageROI.xr || y >= imageROI.yb){ 
          if ((fabs(x)> 2000000000)||(fabs(y)> 2000000000)){
	    cout<<"transform_aligned_track_by_matrix(): pixel outside the image ROI boundaries x:"<<x<<", y:"<<y<<endl;
	    cout<<"transform_aligned_track_by_matrix(): r: "<<r<<endl;
	    cout<<"transform_aligned_track_by_matrix(): transform: "<<transform<<endl;
	  }
	  track[i].image = -1;
	}
	else{
	  track[i].image = image[y*width+x];
	}
      }
    }
  }


  
  void update_aligned_track(vector<AlignedLOLAShot> & track,
			    Eigen::Matrix3d transform, float *image,
			    int width, int height,  struct bbox imageROI)
  {
    transform_aligned_track_by_matrix(track, transform, image, width, height, imageROI);	
    apply_gain_bias(track);
  }
  
  //adjust the position pixel position associated with each point in the track
  //adjut the gain and bias of the reflectance to match the image
  //used!
  void update_aligned_tracks(vector<vector<AlignedLOLAShot> > & tracks,
			     Eigen::Matrix3d transform, float *image,
			     int width, int height, struct bbox imageROI)
  {
    for (unsigned int i = 0; i < tracks.size(); i++){
      transform_aligned_track_by_matrix(tracks[i], transform, image, width, height, imageROI);
    }
    apply_gain_bias(tracks);
  }
 
/*
  vector<Eigen::Matrix3f> load_track_transforms(const std::string& filename)
  {
    vector<Eigen::Matrix3f> l;
    FILE* f = fopen(filename.c_str(), "r");
    
    while (true){
      Eigen::Matrix3f m;
      m(0,0)=1.0;  m(0,1)=0.0;  m(0,2)=0.0;
      m(0,0)=0.0;  m(0,1)=1.0;  m(0,2)=0.0;
      m(0,0)=0.0;  m(0,1)=0.0;  m(0,2)=1.0;
      
      float t[6];
      int ret = fscanf(f, "%g %g %g %g %g %g\n", &t[0], &t[1], &t[2], &t[3], &t[4], &t[5]);
      if (ret != 6){
	return l;
      }
      m(0, 0) = t[0]; m(0, 1) = t[1]; m(0, 2) = t[2];
      m(1, 0) = t[3]; m(1, 1) = t[4]; m(1, 2) = t[5];
      l.push_back(m);
    }
  
    return l;
  }
*/
  void DrawAlignedImageTracks(vector<vector<AlignedLOLAShot> > alignedTracks,
			      GDALDataset* imageDS, string trackType, int thikness)
  {
 
    int imageWidth  = GDALGetRasterXSize(imageDS);
    int imageHeight = GDALGetRasterYSize(imageDS);
    GDALRasterBand* imageBand = imageDS->GetRasterBand(1);
    float *image = new float[imageWidth*imageHeight];
    imageBand->RasterIO(GF_Read, 0, 0, imageWidth, imageHeight,  
			image, imageWidth, imageHeight, 
			GDT_Float32, 0, 0);
    
    float noDataValue = imageBand->GetNoDataValue();
    cout<<"LidarImageAlign.cc:: DrawAlignedImageTracks(): noDataValue="<<noDataValue<<endl;
    
    float minImageVal = 1.0;//imageBand->GetMinimum();
    float maxImageVal = 0.0;//imageBand->GetMaximum();
    for (int row = 0; row < imageHeight; row++){
      for (int col = 0; col < imageWidth; col++){
	if (image[row*imageWidth+col]!=noDataValue){
	  if (image[row*imageWidth+col] > maxImageVal){
	    maxImageVal = image[row*imageWidth+col];
	  }
	  if (image[row*imageWidth+col] < minImageVal){
	    minImageVal = image[row*imageWidth+col];
	  }
        }
      }
    }
    
    double gainDisplay = 1.0/(maxImageVal-minImageVal);
    double biasDisplay = -gainDisplay*minImageVal;
    
    cout<<"LidarImageAlign.cc:: DrawAlignedImageTracks(): minImageVal="<<minImageVal<<", maxImageVal="
	<<maxImageVal<<", gain="<<gainDisplay<<", bias="<<biasDisplay<<", thikness="<<thikness<<endl;
    
    for (int row = 0; row < imageHeight; row++){
      for (int col = 0; col < imageWidth; col++){
	image[row*imageWidth+col] = gainDisplay*image[row*imageWidth+col]+biasDisplay;
      }
    }
    
    //draw the aligned tracks point
    for (unsigned int t = 0; t < alignedTracks.size(); t++){ //for each track
      for (unsigned int s = 0; s < alignedTracks[t].size(); s++){//for each shot	 
        if ((alignedTracks[t][s].valid == 1)){
	  
	  int ptRow = (int)floor(alignedTracks[t][s].image_y);
	  int ptCol = (int)floor(alignedTracks[t][s].image_x);
          int xl = ptCol - thikness; if(xl < 0) {xl = 0;}
	  int xr = ptCol + thikness; if(xr > imageWidth-1) {xr = imageWidth-1;}
	  int yt = ptRow - thikness; if(yt < 0) {yt = 0;}
	  int yb = ptRow + thikness; if(yb > imageHeight-1) {yb = imageHeight-1;}
	  
	  if (trackType.compare("tracks")==0){
	    for (int row = yt; row < yb; row++){
	      for (int col = xl; col < xr; col++){
		image[row*imageWidth+col] = 1.0;
	      }
	    }
	  }
	  if (trackType.compare("reflectance")==0){
	    for (int row = yt; row < yb; row++){
	      for (int col = xl; col < xr; col++){
		if ((col==xl) || (col==xr-1)){
		  image[row*imageWidth+col] = 1.0;
		}
		else{
		  image[row*imageWidth+col] = gainDisplay*alignedTracks[t][s].synth_image+biasDisplay;
		}
	      }
	    }
	  }
	}
      }
    }

    imageBand->RasterIO(GF_Write, 0, 0, imageWidth, imageHeight,  
			image, imageWidth, imageHeight, 
			GDT_Float32, 0, 0);

  }
 void SaveTransformation(Eigen::Matrix3d trans, int numValidPoints, float matchingError, string transfFilename)
  {
    FILE *output = fopen(transfFilename.c_str(), "w");
    if (output != NULL){
      fprintf(output, "%g %g %g %g %g %g\n", trans(0, 0), trans(0, 1), trans(0, 2), trans(1, 0), trans(1, 1), trans(1, 2));
      fprintf(output, "MATCH_ERROR = %g\n", matchingError);
      fprintf(output, "NUM_VALID_POINTS = %d\n", numValidPoints);
      fclose(output);
    }
    else{
      fprintf(stderr, "Failed to open output file %s.\n", transfFilename.c_str());
    }
  }

  void save_track_data( const std::vector<std::vector<AlignedLOLAShot> >& tracks,
			const std::string& filename)
  {
    
    FILE* f = fopen(filename.c_str(), "w");
    if (f == NULL){
      fprintf(stderr, "Could not open track data file %s for writing.\n", filename.c_str());
      return;
    }
    for (unsigned int i = 0; i < tracks.size(); i++){
      int num_valid = 0;
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0){
	  continue;
	}
	num_valid++;
      }
      // ignore small tracks
      if (num_valid < 25){
	continue;
      }
      
      // first print altitude, then distance, then synthetic image, then corresponding image pixel
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0){
	  continue;
	}
	//fprintf(f, "%g", tracks[i][j].LOLAPt[2][2]);
	fprintf(f, "%g", tracks[i][j].LOLAPts[2].z);
	if (j != tracks[i].size() - 1){
	  fprintf(f, " ");
	}
      }
      
      fprintf(f, "\n");
      float last_x = 0.0, last_y = 0.0;
      float total_distance = 0.0;
      
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0){
	  continue;
	}

	//TODO: make sure we use radians and not deg for lon, lat
	Eigen::Vector3f pt;
	pt[0] = tracks[i][j].LOLAPts[2].x;
	pt[1] = tracks[i][j].LOLAPts[2].y;
	pt[2] = tracks[i][j].LOLAPts[2].z;
	Eigen::Vector3f xyz = geographicToCartesian(pt);
	
	if (j != 0){
	  total_distance += sqrt(pow(xyz[0] - last_x, 2) + pow(xyz[1] - last_y, 2));
	}
	
	fprintf(f, "%g", total_distance);
	if (j != tracks[i].size() - 1){
	  fprintf(f, " ");
	}
	last_x = xyz[0];//.x();
	last_y = xyz[1];//.y();
	
      }
      
      fprintf(f, "\n");
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0){
	  continue;
	}
	fprintf(f, "%g", tracks[i][j].synth_image);
	if (j != tracks[i].size() - 1){
	  fprintf(f, " ");
	}
      }
      
      fprintf(f, "\n");
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0){
	  continue;
	}
	fprintf(f, "%g", tracks[i][j].image);
	if (j != tracks[i].size() - 1){
	  fprintf(f, " ");
	}
      }
      fprintf(f, "\n");
    }
    fclose(f);
    
  }

  void SaveReportFile(vector<vector<LOLAShot> > &trackPts, vector<vector<float> >finalTransfArray,  
                        vector<float> errorArray, string reportFilename)
    {
        FILE *fp = fopen(reportFilename.c_str(),"w");
        int numElements = errorArray.size();
        int bestResult = 0;
        float smallestError = 1000000000.0;
	//FIXME!!!
	
        for (int i = 0; i < numElements; i++){
            fprintf(fp,"index=%d d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f e=%f\n", 
                    i,
                    finalTransfArray[i][0], finalTransfArray[i][1], 
                    finalTransfArray[i][2], finalTransfArray[i][3],
                    finalTransfArray[i][4], finalTransfArray[i][5],
                    errorArray[i]);
            if ( (errorArray[i] < smallestError) && (errorArray[i] > 0) ){
                smallestError = errorArray[i]; 
                bestResult = i;
            }    
      
        }
	
        //print the best transform index and smallest error
        fprintf(fp, "bestTransformIndex = %d, smallest error = %f\n", bestResult, errorArray[bestResult]);
        //print the total number of tracks
        fprintf(fp, "numTotalTracks = %d\n", (int)trackPts.size());
    
        //print the number of features per track
        for (unsigned int ti = 0; ti < trackPts.size(); ti++){
            fprintf(fp,"trackID: %d ", ti);
            fprintf(fp,"numShots: %d ", (int)trackPts[ti].size());
            int numFeatures = 0;
            for (unsigned int si = 0; si < trackPts[ti].size(); si++){
                if ((trackPts[ti][si].featurePtRefl == 1.0) && (trackPts[ti][si].valid == 1) && 
                    (trackPts[ti][si].reflectance != 0)) { 
                    numFeatures++;  
                }
            }
            fprintf(fp,"numFeatures: %d\n", numFeatures);
        }


        fclose(fp);
    }


 



 
}
