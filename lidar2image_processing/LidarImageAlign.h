// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef LIDAR_IMAGE_ALIGN_H
#define LIDAR_IMAGE_ALIGN_H

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "../lidar_processing/TracksGCP.h"

namespace at
{
#define DEFAULT_SEARCH_TRANS_WINDOW 20
#define DEFAULT_SEARCH_TRANS_STEP 5.0
#define DEFAULT_SEARCH_THETA_WINDOW (M_PI / 10)
#define DEFAULT_SEARCH_THETA_STEP (M_PI / 40)

#define LIMA 1
#define LIDEM 2

#define MAX_GAUSS_NEWTON_STEPS 100
#define CHUNKSIZE   1

struct CoregistrationParams{
  int matchingMode; //LIMA, LIDEM
  int reflectanceType;//NO, LAMBERT, LUNAR-LAMBERT
  int analyseFlag; 
  int useReflectanceFeatures;
  int topPercentFeatures;
  std::vector<float> samplingStep;
  std::vector<float> matchWindowHalfSize;
  int maxNumIter;
  int maxNumStarts; 
  double noDataVal;
  float minConvThresh;
};

 
/**
 * AlignedLOLAShots are LOLAShots which a transformation can be
 * applied to. The transformed image coordinates and synthetic image value
 * can be stored.
 **/
class AlignedLOLAShot : public LOLAShot
{
public:
  explicit AlignedLOLAShot( LOLAShot& s) : LOLAShot(s), image_x(-1), image_y(-1), synth_image(-1) {};
  ~AlignedLOLAShot(){};

  int   image_x, image_y; // adjusted position in the image
  float synth_image;      // estimated image intensity after applying gain

};
   void DrawAlignedImageTracks(vector<vector<AlignedLOLAShot> > alignedTracks,
			       GDALDataset* imageDS, string trackType, int thikness);
/**
 * Returns a homography from coordinates in aligned to coordinates in aligned2.
 **/
//Matrix3x3 find_track_homography(vector<vector<AlignedLOLAShot> > aligned, vector<vector<AlignedLOLAShot> > aligned2);
  Eigen::Matrix3f find_track_homography(vector<vector<AlignedLOLAShot> > aligned, vector<vector<AlignedLOLAShot> > aligned2);
/**
 * Find a homography from image1 to image2 coordinates. 
 * The variables lonstart, lonend, latstart, latend define the overlapping region between the images.
 **/
//Matrix3x3 find_image_homography(char* image1, char* image2, float lonstart, float lonend, float latstart, float latend);

/**
 * Align a set of LOLA tracks to an image. Returns the transformed points, and the homography found is
 * set as trans. If an image_file is specified, a visualization of the alignment is saved to this file.
 **/
  /*
  vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> > & trackPts,
						  cv::Mat  cubImage, cv::Rect roi,
						  Eigen::Matrix3f & trans);
  */
/**
 * Align a set of LOLA tracks to an image file by searching on an image pyramid.
 * Returns the transformed points, and the homography found is
 * set as trans. If an outputImage is specified, a visualization of the alignment is saved to this file.
 **/
vector<vector<AlignedLOLAShot> > align_to_image_pyramid(vector<vector<LOLAShot> > & trackPts, const string & image_file,
							Eigen::Matrix3d & trans, int *validAlignment, float *error,
							string resultDirname);
/**
 * Compute squared error between a set of tracks synth_image and image.
 * Set numpoints as the number of valid points.
 **/
 float compute_reconstruct_error_ls(vector<AlignedLOLAShot> & track, int* numpoints=NULL);
 float compute_reconstruct_error_ls(vector<vector<AlignedLOLAShot> > & tracks, int* numpoints=NULL);
 Eigen::MatrixXf compute_reconstruct_error_vec(vector<AlignedLOLAShot> & track, int num_points);

/**
 * Find translation and rotation to align tracks to image using brute force search.
 * The translation and rotation to search over can be specified.
 * Returns the translation and rotation as an affine transform.
 **/
//Matrix3x3 find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, ImageView<PixelGray<float> > & cub,
//		int transSearchWindow=DEFAULT_SEARCH_TRANS_WINDOW, int transSearchStep=DEFAULT_SEARCH_TRANS_STEP, 
//		float thetaSearchWindow=DEFAULT_SEARCH_THETA_WINDOW, float thetaSearchStep=DEFAULT_SEARCH_THETA_STEP);
/*
 Eigen::Matrix3f find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, cv::Mat cubImage,
				       int transSearchWindow=DEFAULT_SEARCH_TRANS_WINDOW, 
				       int transSearchStep=DEFAULT_SEARCH_TRANS_STEP, 
				       float thetaSearchWindow=DEFAULT_SEARCH_THETA_WINDOW, 
				       float thetaSearchStep=DEFAULT_SEARCH_THETA_STEP);
*/
/**
 * Find a homography to align tracks to image.
 * Can specify a starting point in the search with matrix.
 **/
//Matrix3x3 gauss_newton_homography(vector<AlignedLOLAShot> & track, ImageView<PixelGray<float> > cubImage,
//		Matrix3x3 matrix=Matrix3x3(1,0,0,0,1,0,0,0,1));
 //Eigen::Matrix3f gauss_newton_homography(vector<AlignedLOLAShot> & track, cv::Mat cubImage,
 //					 Eigen::Matrix3f matrix/*=Eigen::Matrix3f(1,0,0,0,1,0,0,0,1)*/);
/**
 * Find an affine transform to align tracks to image. Runs faster than finding a homography.
 * Can specify a starting point in the search with matrix.
 **/
//Matrix3x3 gauss_newton_affine(vector<AlignedLOLAShot> & track, ImageView<PixelGray<float> > cubImage,
//		Matrix3x3 matrix=Matrix3x3(1,0,0,0,1,0,0,0,1));
 //Eigen::Matrix3f gauss_newton_affine(vector<AlignedLOLAShot> & track,
 //				     cv::Mat cubImage, cv::Rect roi,
 //				     Eigen::Matrix3f matrix);
 Eigen::Matrix3f gauss_newton_affine(vector<AlignedLOLAShot> & track,
				     float* image, int width, int height,
				     struct bbox imageROI, Eigen::Matrix3f matrix, float *r_err);
 Eigen::Matrix3d gaussNewtonAffine(vector<AlignedLOLAShot> & track,
				     float* image, int width, int height,
				     struct bbox imageROI, Eigen::Matrix3d matrix, float *r_err);
 /**
 * Compute gain and bias from a vector of LOLAShots. These factors are
 * applied to the reflectance to create a synthetic image. We choose
 * values such that the sum of the synthetic pixels equals the sum of
 * the original pixels. Currently the bias is 0.
 **/
 Eigen::Vector2f ComputeGainBiasFactor( const std::vector<AlignedLOLAShot>& );
 Eigen::Vector2f ComputeGainBiasFactor( const std::vector<std::vector<AlignedLOLAShot> >& );
/**
 * Transform track coordinates and image samples based on the matrix transform.
 *
 * Also recomputes synthetic image by computing new gain and bias from image.
 **/
//void transform_track(vector<AlignedLOLAShot> & track, Matrix3x3 transform, ImageView<PixelGray<float> >& cub);
//void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, ImageView<PixelGray<float> >& cub);
//void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, string cubFile);

 //void update_track(vector<AlignedLOLAShot> & track, Eigen::Matrix3f transform, cv::Mat cub, cv::Rect roi);
 //void update_tracks(vector<vector<AlignedLOLAShot> > & tracks, Eigen::Matrix3f transform, cv::Mat cub, cv::Rect roi);
 void update_aligned_track(vector<AlignedLOLAShot> & track,Eigen::Matrix3d transform, float *image,
			   int width, int height,  struct bbox imageROI);
 void update_aligned_tracks(vector<vector<AlignedLOLAShot> > & tracks, Eigen::Matrix3d transform, float *image,
			    int width, int height, struct bbox imageROI);
 //void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Eigen::Matrix3f transform, string cubFilename);
/**
 * Transform underlying coordinates of tracks. New calls to transform track coordinates
 * will apply their transform to the coordinates resulting from the call to these functions
 * rather than to the original coordintaes.
 *
 * Can specify a single matrix to transform all tracks by, or per-track transformations.
 **/
//void transform_tracks_by_matrices(vector<vector<LOLAShot> > & tracks, vector<Matrix3x3> matrices);
//void transform_tracks_by_matrix(vector<vector<LOLAShot> > & tracks, Matrix3x3 M);
 void transform_tracks_by_matrices(vector<vector<LOLAShot> > & tracks, vector<Eigen::Matrix3d> matrices);
 void transform_tracks_by_matrix(vector<vector<LOLAShot> > & tracks, Eigen::Matrix3d M,  bbox imageROI);
 // void transform_track_by_matrix(vector<AlignedLOLAShot> & track, Eigen::Matrix3f transform, cv::Mat image, cv::Rect roi);
/**
 * Save track data in format suitable for visualization.
 **/
void save_track_data( const std::vector<std::vector<AlignedLOLAShot> >&, const std::string& filename);
 void SaveTransformation(Eigen::Matrix3d trans, int validAlignment, float matchingError, std::string transFilename);
/**
 * Load list of affine transformation matrices from file.
 * Each line in the file has the six entries of one matrix separated by spaces.
 **/
//vector<Matrix3x3> load_track_transforms(const std::string& filename);
 vector<Eigen::Matrix3f> load_track_transforms(const std::string& filename);
/**
 * Create vector of AlignedLOLAShots (used for transforms) from vector of LOLAShots.
 **/
std::vector<std::vector< AlignedLOLAShot> > initialize_aligned_lola_shots(std::vector<std::vector<LOLAShot> >& tracks);
/*
* Draw aligned tracks on image. The tracks are shown in their synthetic image intensity
*/
//cv::Mat DrawAlignedImageTracks( vector<vector<AlignedLOLAShot> >, cv::Mat image,
//				 float gain, float bias, int thikness);
// cv::Mat DrawAlignedTracks( vector<vector<AlignedLOLAShot> > alignedTracks, cv::Mat image,
//			    float gain, float bias, int thikness);

 void UpdateGCP( const std::vector<std::vector<AlignedLOLAShot> >& trackPts, 
		 const std::string& camCubFile, 
		 vector<gcp>&  gcpArray);


}
#endif//LIDAR_IMAGE_ALIGN_H
