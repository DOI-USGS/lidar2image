/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

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
  
#if 0
  //computes the least square error between the image and the synthetic image
  float compute_reconstruct_error_ls(vector<AlignedLOLAShot> & track, int* numPoints)
  {
    float err = 0.0;
    int num_points = 0;
    
    for (unsigned int j = 0; j < track.size(); j++){
      if ((track[j].image != -1) && (track[j].valid == 1)){
	err += pow(track[j].synth_image - track[j].image, 2);
	num_points++;
      }
    }
    
    cout<<"compute_reconstruct_error_ls (single track): num_points="<<num_points<<endl;

    if (numPoints != NULL){
      *numPoints = num_points;
    }
    if (num_points == 0){
      return -1.0;
    }

    cout<<"compute_reconstruct_error_ls (single track): num_points="<<num_points<<endl;
    return err / num_points;
  }
  
#endif
#if 0
  float compute_reconstruct_error_ls(vector<vector<AlignedLOLAShot> > & tracks, int* numPoints)
  {
    float error = 0.0;
    int numValidPoints = 0;

    cout<<"compute_reconstruct_error_ls (multi-track)"<<endl;
    
    for (unsigned int i = 0; i < tracks.size(); i++){
      /*
      float numValidPointsInTrack;
      float normErrorInTrack = compute_reconstruct_error_ls(tracks[i], &numValidPointsInTrack);
      
      //update the number of valid points
      numValidPoints = numValidPoints + numValidPointsInTrack;
      
      //update the error
      float errorInTrack = 0;
      if (normErrorInTrack!=-1){
	errorInTrack = normErrorInTrack*numValidPointsInTrack;
      } 
      error = error + errorInTrack;
    }
    */
    
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if ((tracks[i][j].image != -1) && (tracks[i][j].valid == 1)){	  
	  error += pow(tracks[i][j].synth_image - tracks[i][j].image, 2);
	  numValidPoints++;
	}
      }
    }

    if (numPoints != NULL){
      *numPoints = numValidPoints;
    }
    if (numValidPoints == 0){
      return -1;
    }
    
    return error / numValidPoints;
  }
#endif
/*
  float compute_reconstruct_error_ls(vector<vector<AlignedLOLAShot> > & tracks, int* numPoints)
  {
    float err = 0.0;
    int num_points = 0;

    cout<<"compute_reconstruct_error_ls (multi-track)"<<endl;

    for (unsigned int i = 0; i < tracks.size(); i++){
      for (unsigned int j = 0; j < tracks[i].size(); j++){
	if ((tracks[i][j].image != -1) && (tracks[i][j].valid == 1)){	  
	  err += pow(tracks[i][j].synth_image - tracks[i][j].image, 2);
	  num_points++;
	}
      }
    }
    
    if (numPoints != NULL){
      *numPoints = num_points;
    }
    if (num_points == 0){
      return -1;
    }
    
    return err / num_points;
  }
*/
#if 0
  //returns a vector of size numPoints. each element of this vector is
  //the relative error between the observed and reconstructed image intensity
  Eigen::MatrixXf compute_reconstruct_error_vec(vector<AlignedLOLAShot> & track, int numPoints)
  {
    Eigen::MatrixXf r(numPoints, 1);
    int validPtIndex = 0;


    for (unsigned int i = 0; i < track.size(); i++){

      //here is a problem. what happens if the condition below is not satisfied.
      cout<<"validPtIndex="<<validPtIndex<<", numPoints="<<numPoints<<endl; 
      /*
      if (validPtIndex < numPoints){
	r(validPtIndex, 0) = 10.0;
      }
      */
      if ((track[i].image != -1) && (track[i].valid == 1) && (validPtIndex < numPoints)){	
	r(validPtIndex, 0) = track[i].image - track[i].synth_image;
	if (fabs(r(validPtIndex, 0)) > 1.0){
	  r(validPtIndex, 0) = 1.0;
	}
        validPtIndex++;    
      }
      
      if (validPtIndex > numPoints){
	cout<<"LidarImageAlign.cc:: compute_reconstruct_error_vec(): ERROR"<<endl;
	cout<<"i="<<i<<", validPtIndex="<<validPtIndex<<" >= numPoints="<<numPoints
            <<", totalNumPoints="<<track.size()<<endl;
      }
     
    }
    cout<<"LidarImageAlign.cc:: compute_reconstruct_error_vec(): numValidPoints="<<validPtIndex
        <<", totalNumPoints="<<track.size()<<endl;
    return r;
  }
#endif

 
  /*
    // brute force search
  Matrix3x3 find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, cv::Mat cub,
                                    int transSearchWindow, int transSearchStep, float thetaSearchWindow, float thetaSearchStep)
    {
        float mid_x = cub.cols() / 2;
        float mid_y = cub.rows() / 2;
        Matrix3x3 center(1, 0, -mid_x, 0, 1, -mid_y, 0, 0, 1);
        Matrix3x3 uncenter(1, 0, mid_x, 0, 1, mid_y, 0, 0, 1);

        Matrix3x3 best;
        float best_score = INFINITY;


        for (float tt = -thetaSearchWindow; tt <= thetaSearchWindow; tt += thetaSearchStep)
        {
            Matrix3x3 purerot(cos(tt), -sin(tt), 0, sin(tt), cos(tt), 0, 0, 0, 1);
            Matrix3x3 rot = uncenter * purerot * center;
            for (int xt = -transSearchWindow; xt <= transSearchWindow; xt += transSearchStep)
                for (int yt = -transSearchWindow; yt <= transSearchWindow; yt += transSearchStep)
                {
                    Matrix3x3 trans(1, 0, xt, 0, 1, yt, 0, 0, 1);
                    transform_tracks(tracks, trans * rot, cub);
                    float score = compute_transform_error(tracks);
                    if (score < best_score)
                    {
                        //printf("%g %d %d %g\n", tt, xt, yt, score);
                        best_score = score;
                        best = trans * rot;
                    }
                }
            //printf("%g\n", tt);
        }
        //printf("Best x: %g y: %g theta: %g\n", best_trans_x, best_trans_y, best_rot);
	
        transform_tracks(tracks, best, cub);

        return best;
    }
  */

#if 0
  // brute force search
  Eigen::Matrix3f find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, cv::Mat cub,
					int transSearchWindow, int transSearchStep,
					float thetaSearchWindow, float thetaSearchStep)
  {
    float mid_x = cub.cols / 2;
    float mid_y = cub.rows / 2;

    Eigen::Matrix3f center;
    center(0,0) = 1; center(0,1) = 0; center(0,2) = -mid_x;
    center(1,0) = 0; center(1,1) = 1; center(1,2) = -mid_y;
    center(2,0) = 0; center(2,1) = 0; center(2,2) = 1;
    Eigen::Matrix3f uncenter;
    uncenter(0,0) = 1; uncenter(0,1) = 0; uncenter(0,2) = mid_x; 
    uncenter(1,0) = 0; uncenter(1,1) = 1; uncenter(1,2) = mid_y; 
    uncenter(2,0) = 0; uncenter(2,1) = 0; uncenter(2,2) = 1;
    
    Eigen::Matrix3f best;
    float best_score = INFINITY;
    cv::Rect roi;
    
    for (float tt = -thetaSearchWindow; tt <= thetaSearchWindow; tt += thetaSearchStep){
      //Matrix3x3 purerot(cos(tt), -sin(tt), 0, sin(tt), cos(tt), 0, 0, 0, 1);
      Eigen::Matrix3f purerot;
      purerot(0,0) = cos(tt);  purerot(0,1) =-sin(tt); purerot(0,2) = 0; 
      purerot(1,0) = sin(tt);  purerot(1,1) = cos(tt); purerot(1,2) = 0; 
      purerot(2,0) =   0;      purerot(2,1) = 0;       purerot(2,2) = 1;
      
      Eigen::Matrix3f rot = uncenter * purerot * center;
 
      for (int xt = -transSearchWindow; xt <= transSearchWindow; xt += transSearchStep){
	for (int yt = -transSearchWindow; yt <= transSearchWindow; yt += transSearchStep){
	  
	  Eigen::Matrix3f trans;
	  trans(0,0) = 1; trans(0,1) = 0; trans(0,2) = xt; 
	  trans(1,0) = 0; trans(1,1) = 1; trans(1,2) = yt; 
	  trans(1,1) = 0; trans(2,1) = 0; trans(2,2) = 1;

	  update_tracks(tracks, trans * rot, cub, roi);

	  float score = compute_reconstruct_error_ls(tracks);

	  if (score < best_score){
	    //printf("%g %d %d %g\n", tt, xt, yt, score);
	    best_score = score;
	    best = trans * rot;
	  }
	}
	//printf("%g\n", tt);
      }
    }
    //printf("Best x: %g y: %g theta: %g\n", best_trans_x, best_trans_y, best_rot);
    //cv::Rect roi; 
    update_tracks(tracks, best, cub, roi);
    
    return best;
  }

  #endif
#if 0
  Eigen::MatrixXf compute_jacobian_affine(vector<AlignedLOLAShot> & track,
					  float* image, int width, int height, 
                                          int numPoints)
  {
    Eigen::MatrixXf J(numPoints, 6); 
    int index = 0;

    for (unsigned int i = 0; i < track.size(); i++){
      // TODO: we may ignore some points that weren't in the image when num_points was computed,
      // but are with the new transformation
      if (track[i].image == -1 || track[i].valid==0 || index >= numPoints){
	continue;
      }
      
      int x = track[i].image_x;
      int y = track[i].image_y;
      
      double h = 1.0;
      
      //fininte difference with three points to right and left
      double dx, dy;
      if (x >= width - 3 || x <= 3 || y >= height - 3 || y <= 3){
	dx = 0.0; dy = 0.0;
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
    
    return J;
  }

#endif

 
#if 0
  //Gauss Newton solution
  //returns the matching error r_err; 
  Eigen::Matrix3f gauss_newton_affine(vector<AlignedLOLAShot> & track,
				      float* image, int width, int height,
				      struct bbox imageROI, Eigen::Matrix3f matrix,                                      
				      float *r_err)
  {
    Eigen::Matrix3f B = matrix;
    
    int numInitValidPoints = 0;
    float curr_err = 100000.0;
    float prev_err = compute_reconstruct_error_ls(track, &numInitValidPoints);
    cout<<"LidarImageAlign.cc:: gauss_newton_affine(): numInitValidPoints="<<numInitValidPoints<<endl;
    
    if (numInitValidPoints <= 0){ // no points
      *r_err = curr_err;
      return matrix;
    }

    cout<<"LidarImageAlign.cc:: gauss_newton_affine(): B="<<endl;
    cout<<B<<endl;
    
    for (int count = 0; count < MAX_GAUSS_NEWTON_STEPS; count++){
     
      Eigen::MatrixXf J  = compute_jacobian_affine(track, image, width,
						   height, numInitValidPoints);

      Eigen::MatrixXf JT = J.transpose();
      Eigen::MatrixXf r  = compute_reconstruct_error_vec(track, numInitValidPoints);

      cout<<"LidarImageAlign.cc:: gauss_newton_affine(): sqn of r = "<<r.squaredNorm()<<endl;

      Eigen::JacobiSVD<Eigen::MatrixXf> svd(JT*J, Eigen::ComputeThinU | Eigen::ComputeThinV);
      const Eigen::VectorXf s  = svd.singularValues();
      const Eigen::MatrixXf U  = svd.matrixU();
      const Eigen::MatrixXf V  = svd.matrixV();

      //cout<<"s="<<s<<endl;
      
      Eigen::MatrixXf S(6, 6);
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

      cout<<"LidarImageAlign.cc:: gauss_newton_affine(): S="<<endl;
      cout<<S<<endl;

      Eigen::MatrixXf pseudoinverse = V * S * U.transpose();
      
      Eigen::MatrixXf delta = - pseudoinverse * (JT * r);
      
      cout<<"LidarImageAlign.cc:: gauss_newton_affine(): delta="<<endl;
      cout<<delta<<endl;

      Eigen::Matrix3f last_matrix = B;

      B(0, 0) += delta(0, 0); B(0, 1) += delta(1, 0); B(0, 2) += delta(2, 0);
      B(1, 0) += delta(3, 0); B(1, 1) += delta(4, 0); B(1, 2) += delta(5, 0);

      cout<<"LidarImageAlign.cc:: gauss_newton_affine(): B="<<endl;
      cout<<B<<endl;
      
      update_aligned_track(track, B, image, width, height, imageROI);
      
      curr_err = compute_reconstruct_error_ls(track);

      cout<<"LidarImageAlign.cc:: gauss_newton_affine: err="<<prev_err<<" before iter "<<count<<endl;
      cout<<"LidarImageAlign.cc:: gauss_newton_affine: err ="<<curr_err<<" after iter  "<<count<<endl;

      if (curr_err == -1){
	fprintf(stderr, "Matrix jumped too far, giving up.\n");
        *r_err = prev_err;  
	return last_matrix;
      }

      if (curr_err > prev_err){
	B = last_matrix;
        curr_err = prev_err; 
	break;
      }
      prev_err = curr_err;
    }

    *r_err = curr_err;
   
     return B;

  }
#endif
 
#if 0
  //not used
  Matrix3x3 find_image_homography(char* image1, char* image2, float lonstart, float lonend, float latstart, float latend)
    {
        camera::IsisCameraModel model1(image1);
        camera::IsisCameraModel model2(image2);

        vector<vector<AlignedLOLAShot> > aligneda, alignedb;
        aligneda.push_back(vector<AlignedLOLAShot>());
        alignedb.push_back(vector<AlignedLOLAShot>());
        for (float x = lonstart; x <= lonend; x += 0.1)
            for (float y = latstart; y <= latend; y += 0.1)
            {
                // height at one point, we treat the surface as a plane
	      Vector3 lon_lat_rad( x, y, 1735.0 * 1000);
	      Vector3 xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
	      Vector2 a = model1.point_to_pixel(xyz);
	      Vector2 b = model2.point_to_pixel(xyz);

	      std::vector<pointCloud> p;
	      LOLAShot l1(p), l2(p);
	      AlignedLOLAShot s1(l1), s2(l2);
	      s1.image_x = a.x(); s1.image_y = a.y();
	      s2.image_x = b.x(); s2.image_y = b.y();
	      aligneda[0].push_back(s1);
	      alignedb[0].push_back(s2);
            }
        Matrix3x3 H = find_track_homography(aligneda, alignedb);
        return H;
    }
#endif
 // returns a transformation matrix from coordinates in image 2 to coordinates in image 1
  Eigen::Matrix3f find_track_homography(vector<vector<AlignedLOLAShot> > aligned, vector<vector<AlignedLOLAShot> > aligned2)
  {
    int num_points = 0;
    for (unsigned int i = 0; i < aligned.size(); i++)
      {
	for (unsigned int j = 0; j < aligned[i].size(); j++)
	  if (aligned[i][j].image_x >= 0 && aligned2[i][j].image_x >= 0)
	    num_points++;
      }
    // [a b c; d e f; 0 0 1] * [x;y;1] = [x';y';1] = [ax + by + c; dx + ey + f; 1]
    // [x y 1 0 0 0; 0 0 0 x y 1] * [a;b;c;d;e;f]
    // [a b c; d e f; g h i] * [x;y;1] = [x';y';w'] = [ax + by + c; dx + ey + f; gx + hy + i]
    // [x y 1 0 0 0 0 0 0; 0 0 0 x y 1 0 0 0; 0 0 0 0 0 0 x y 1] * [a;b;c;d;e;f;g;h;i]
    // [h11 x + h12 y + h13 - 
    // Matrix<double> A(2 * num_points, 9);
    // Matrix<double> b(2 * num_points, 1);

    Eigen::MatrixXf A(2*num_points, 9);
    Eigen::MatrixXf b(2*num_points, 1);

    int index = 0;
    for (unsigned int i = 0; i < aligned.size(); i++)
      for (unsigned int j = 0; j < aligned[i].size(); j++)
	if (aligned[i][j].image_x >= 0 && aligned2[i][j].image_x >= 0)
	  {
	    A(2 * index    , 0) = aligned2[i][j].image_x;
	    A(2 * index    , 1) = aligned2[i][j].image_y;
	    A(2 * index    , 2) = 1;
	    A(2 * index    , 3) = 0;
	    A(2 * index    , 4) = 0;
	    A(2 * index    , 5) = 0;
	    A(2 * index    , 6) = -aligned[i][j].image_x * aligned2[i][j].image_x;
	    A(2 * index    , 7) = -aligned[i][j].image_x * aligned2[i][j].image_y;
	    A(2 * index    , 8) = -aligned[i][j].image_x;
	    A(2 * index + 1, 0) = 0;
	    A(2 * index + 1, 1) = 0;
	    A(2 * index + 1, 2) = 0;
	    A(2 * index + 1, 3) = aligned2[i][j].image_x;
	    A(2 * index + 1, 4) = aligned2[i][j].image_y;
	    A(2 * index + 1, 5) = 1;
	    A(2 * index + 1, 6) = -aligned[i][j].image_y * aligned2[i][j].image_x;
	    A(2 * index + 1, 7) = -aligned[i][j].image_y * aligned2[i][j].image_y;
	    A(2 * index + 1, 8) = -aligned[i][j].image_y;
	    b(2 * index    , 0) = 0;
	    b(2 * index + 1, 0) = 0;
	    //printf("(%d %d) (%d %d)\n", aligned2[i][j].image_x, aligned2[i][j].image_y, aligned[i][j].image_x, aligned[i][j].image_y);
	    index++;
	  }
   
    Eigen::MatrixXf U (2*num_points, 9);
    Eigen::MatrixXf S (9,9);
    Eigen::MatrixXf VT(9,9);
    Eigen::VectorXf s(9); 

    //Matrix<double> U(2 * num_points, 9), S(9, 9), VT(9, 9);
    //Vector<double> s(9);

    //FIXME!!!!
    // svd(A, U, s, VT);

    // Matrix3x3 H(VT(8, 0), VT(8, 1), VT(8, 2), VT(8, 3), VT(8, 4), VT(8, 5), VT(8, 6), VT(8, 7), VT(8, 8));
    
    Eigen::Matrix3f H;
    H(0,0)=VT(8,0); H(0,1)=VT(8,1); H(0,2)=VT(8,2);
    H(1,0)=VT(8,3); H(1,1)=VT(8,4); H(1,2)=VT(8,5);
    H(2,0)=VT(8,6); H(2,1)=VT(8,7); H(2,2)=VT(8,8);

    return H;
  }

  #if 0
     // J_f(B)_{ij} = dF_i / dB_j
  // transform may be marked as invalid if points have left the image window
  Eigen::MatrixXf compute_jacobian_homography(vector<AlignedLOLAShot> & track, Eigen::Matrix3f B, cv::Mat img, int num_points)
    {
      /*
    	InterpolationView<EdgeExtensionView<ImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation> imgp
    		= interpolate(img, BilinearInterpolation(), ConstantEdgeExtension());
      */
        //Matrix<double> J(num_points, 9);
	Eigen::MatrixXf J(num_points, 9);
        int index = 0;
        for (unsigned int i = 0; i < track.size(); i++)
        {
            // TODO: we may ignore some points that weren't in the image when num_points was computed,
            // but are with the new transformation
            if (track[i].image == -1 || track[i].synth_image == -1 || index >= num_points)
                continue;

            double ox = (double)track[i].imgPt[2].x, oy = (double)track[i].imgPt[2].y;
            int x = track[i].image_x;
            int y = track[i].image_y;

            double h = 1.0;
		
            // fininte difference with three points to right and left
            double dx=0.0, dy =0.0, dz = 0.0, h31 = 1.0, h32 = 1.0, h33 = 1.0;

	    /*	    
            if (x >= img.rows() - 3 || x <= 3 || y >= img.cols() - 3 || y <= 3)
                dx = 0.0, dy = 0.0, dz = 0.0, h31 = 1.0, h32 = 1.0, h33 = 1.0;
            else
            {
                dx = (-1.0/60) * img(x-3, y) + (3.0/20) * img(x-2, y) - (3.0 / 4) * img(x-1, y) +
                    (3.0/4) * img(x+1, y) - (3.0 / 20) * img(x+2, y) + (1.0 / 60) * img(x+3, y);
                dy = (-1.0/60) * img(x, y-3) + (3.0/20) * img(x, y-2) - (3.0 / 4) * img(x, y-1) +
                    (3.0/4) * img(x, y+1) - (3.0 / 20) * img(x, y+2) + (1.0 / 60) * img(x, y+3);
                h = B(2, 0) * ox + B(2, 1) * oy + B(2, 2);
                if (ox == 0 || oy == 0 || x == 0 || y == 0)
                    continue;
                float m = (float)y / x;
                if (x > y)
                {
                    float a = (B(0, 0) * ox + B(0, 1) * oy + B(0, 2)) / (1 + (double)x);
                    dz = (-1.0/60) * imgp(x-3, y - 3*m) + (3.0/20) * imgp(x-2, y-2*m) - (3.0 / 4) * imgp(x-1, y-m) +
                        (3.0/4) * imgp(x+1, y+m) - (3.0 / 20) * imgp(x+2, y+2*m) + (1.0 / 60) * imgp(x+3, y+3*m);
                    h31 = 1 / ox * (a - B(2, 1) * oy - B(2, 2)) - B(2, 0);
                    h32 = 1 / oy * (a - B(2, 0) * ox - B(2, 2)) - B(2, 1);
                    h33 = (a - B(2, 0) * ox - B(2, 1) * oy) - B(2, 2);
                }
                else
                {
                    float a = (B(1, 0) * ox + B(1, 1) * oy + B(1, 2)) / (1 + (double)y);
                    dz = (-1.0/60) * imgp(x-3/m, y - 3) + (3.0/20) * imgp(x-2/m, y-2) - (3.0 / 4) * imgp(x-1/m, y-1) +
                        (3.0/4) * imgp(x+1/m, y+1) - (3.0 / 20) * imgp(x+2/m, y+2) + (1.0 / 60) * imgp(x+3/m, y+3);
                    h31 = 1 / ox * (a - B(2, 1) * oy - B(2, 2)) - B(2, 0);
                    h32 = 1 / oy * (a - B(2, 0) * ox - B(2, 2)) - B(2, 1);
                    h33 = (a - B(2, 0) * ox - B(2, 1) * oy) - B(2, 2);
                }
	    
                if (h31 == 0.0 || h32 == 0.0 || h33 == 0.0)
                {
                    h31 = 1.0; h32 = 1.0; h33 = 1.0; dz = 0.0; // don't divide by zero
                }
            }
	    */
            // compute the derivatives
            J(index, 0) = dx / (h / track[i].imgPt[2].x);
            J(index, 1) = dx / (h / track[i].imgPt[2].y);
            J(index, 2) = dx / h;
            J(index, 3) = dy / (h / track[i].imgPt[2].x);
            J(index, 4) = dy / (h / track[i].imgPt[2].y);
            J(index, 5) = dy / h;
            J(index, 6) = dz / h31;
            J(index, 7) = dz / h32;
            J(index, 8) = dz / h33;
            index++;

        }
	    
        return J;
    }
#endif
 #if 0 
     /**
     * B = transformation matrix
     * f_i(B) = synthetic image at track point i given B
     * r_i(B) = error at track point i given B
     * B_{s+1} = B_s + \delta
     * (J_f^T J_f) \delta = J_f^T r
     * Must compute J_f(B)
     * J_f(B)_{ij} = dF_i / dB_j
     **/
    // assumes we have already called transform_track to set image and synth_image
  Eigen::Matrix3f gauss_newton_homography(vector<AlignedLOLAShot> & track, cv::Mat cubImage,
					  Eigen::Matrix3f matrix)
  {
    
    //Matrix3x3 B = matrix;
    Eigen::Matrix3f B = matrix;
    
    int num_points = 0;
    float err, last_err = compute_reconstruct_error_ls(track, &num_points);
    if (num_points <= 0) // no points
      return matrix;

    cv::Rect roi;
    for (int count = 0; count < MAX_GAUSS_NEWTON_STEPS; count++){

	Eigen::MatrixXf J = compute_jacobian_homography(track, B, cubImage, num_points);
	Eigen::MatrixXf trans = J.transpose();
	Eigen::MatrixXf  r = compute_reconstruct_error_vec(track, num_points);
	
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(trans*J, Eigen::ComputeThinU | Eigen::ComputeThinV);
	const Eigen::VectorXf s = svd.singularValues();
	const Eigen::MatrixXf U = svd.matrixU();
	const Eigen::MatrixXf V = svd.matrixV();

	Eigen::MatrixXf S(6, 6);
	for (int i = 0; i < 9; i++){
	    if (s(i) > 10e-15) // don't divide by 0
	      S(i, i) = 1.0 / s(i);
	}
	
	//Matrix<double> pseudoinverse = transpose(VT) * S * transpose(U);
	//Matrix<double> delta = pseudoinverse * -(trans * r);
	Eigen::MatrixXf pseudoinverse;
        Eigen::MatrixXf delta;
        Eigen::Matrix3f last_matrix = B;
	B(0, 0) += delta(0, 0);
	B(0, 1) += delta(1, 0);
	B(0, 2) += delta(2, 0);
	B(1, 0) += delta(3, 0);
	B(1, 1) += delta(4, 0);
	B(1, 2) += delta(5, 0);
	B(2, 0) += delta(6, 0);
	B(2, 1) += delta(7, 0);
	B(2, 2) += delta(8, 0);
	update_track(track, B, cubImage, roi);
	err = compute_reconstruct_error_ls(track);
	if (err == -1)
	  {
	    fprintf(stderr, "Matrix jumped too far, giving up.\n");
	    return last_matrix;
	  }
	//printf("Error: %g Old Error: %g\n", err, last_err);
	if (err >= last_err)
	  {
	    B = last_matrix;
	    break;
	  }
	last_err = err;
      }
    
    return B;
  }
#endif


 
}
