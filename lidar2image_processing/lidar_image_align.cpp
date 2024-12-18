/**
 * This program takes as input a cub file which is an image of a planetary surface,
 * and a CSV file of LOLA lidar points. It performs a brute force search over a
 * local area to find the transform which best fits the lidar points to the image
 * according to their estimate reflectance.
 **/

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>

#include "../common/StringUtils.h"
#include "../lidar_processing/TracksLOLA.h"
#include "../lidar_processing/TracksFeatures.h"
#include "LidarImageAlign.h"


using namespace std;
using namespace at;

int main( int argc, char *argv[] )
{
  // command line options
  string LOLAFilename;
  string imageFilename;
  string settingsFilename;
  string resultDirname="";

  int LOLAFilenameFound = 0;
  int imageFilenameFound = 0;
  int settingsFilenameFound = 0;
  int resultDirnameFound = 0;
  int testFlag = 0;
  int testError = 0;
  
  for (int i = 1; i < argc; i++) { //We will iterate over argv[] to get the parameters stored inside.
   //Note that we're starting on 1 because we don't need to know the
   //path of the program, which is stored in argv[0]
   if (i + 1 != argc) { // Check that we haven't finished parsing already

     if (string(argv[i]).compare("-t")==0) {
       testFlag = 1;
     }
     else if (string(argv[i]).compare("-l")==0) {
       LOLAFilename = argv[i + 1];
       LOLAFilenameFound = 1;
     }
     else if (string(argv[i]).compare("-i")==0) {
       imageFilename = argv[i + 1];
       imageFilenameFound = 1;
     }
     else if (string(argv[i]).compare("-s")==0) {
       settingsFilename = argv[i + 1];
       settingsFilenameFound = 1; 
     }
     else if (string(argv[i]).compare("-r")==0) {
       resultDirname = argv[i + 1];
       resultDirnameFound = 1; 
     }
   }
 }

 int usageError = 0;
 
 if (LOLAFilenameFound == 0){
   cout<<"Error! LOLA filename not found"<<endl;
   usageError = 1;
 }
 if (imageFilenameFound==0){
   cout<<"Error! Image filename or image pyramid not found"<<endl;
   usageError = 1;
 }

 if (usageError == 1){
   std::cout << "Usage is "<< argv[0]<< " -l LOLAFilename -i imageFilename "
	     <<" [-t] [-s settingsFilename] [-r resultDirname]"<<endl;
   std::cin.get();
   exit(0);
 }

 vector<vector<LOLAShot> > trackPts;
 trackPts = LOLAFilesRead(LOLAFilename);

 Eigen::Matrix3d trans;
 trans(0,0)=1; trans(0,1)=0; trans(0,2)=0;
 trans(1,0)=0; trans(1,1)=1; trans(1,2)=0;
 trans(2,0)=0; trans(2,1)=0; trans(2,2)=1;
 
 // use this later, only if saving to file
 vector<gcp> gcpArray;
 gcpArray = ComputeSalientLOLAFeatures(trackPts);

 int validAlignment = 1;
 float matchingError = 0.0;
 vector<vector< AlignedLOLAShot> > aligned;
 aligned = align_to_image_pyramid(trackPts, imageFilename, trans, &validAlignment, &matchingError, resultDirname);

 if (resultDirname.compare("")!=0){
   
   string gcpFilename;
   string transfFilename;
   string alignedTracksFilename;
    
   transfFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))+"_transf.txt";
   SaveTransformation(trans, validAlignment, matchingError, transfFilename);
 
   UpdateGCP(aligned, imageFilename, gcpArray);
   //gcpFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename));
   //SaveGCPoints(gcpArray, gcpFilename);
    
   gcpFilename = resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))+"_gcp.txt";
   cout<<"gcpFilename="<<gcpFilename<<endl;
   SaveGCP(gcpArray, gcpFilename);
   
   alignedTracksFilename=resultDirname+"/"+GetFilenameNoPath(GetFilenameNoExt(imageFilename))+"_aligned_tracks.txt";
   save_track_data(aligned, alignedTracksFilename);
 }
 
 if (testFlag==1){
   
   Eigen::Matrix3d transTest;
   string alignedTracksFilename=GetFilenameNoExt(imageFilename)+"_transf.txt";
   ifstream file (alignedTracksFilename.c_str());
   std::string line;
   if (file.is_open()){
     std::getline(file, line);
     //cout<<line<<endl;
     stringstream lineStream(line);
     lineStream>>transTest(0,0)>>transTest(0,1)>>transTest(0,2)>>transTest(1,0)>>transTest(1,1)>>transTest(1,2);  
     transTest(2,0) = 0; transTest(2,1) = 0; transTest(2,2) = 1;                
     file.close();
   }
   else{
     cout<<"file not found "<< alignedTracksFilename<<endl;
   }

   cout<<"lidar_image_align.cc:: main(): transTest="<<endl;
   cout<<transTest<<endl;
   
   double distance = 0;
   for (int row = 0; row<2; row++){
     for (int col = 0; col<3; col++){
       distance = distance + (trans(row,col)-transTest(row, col))* (trans(row,col)-transTest(row, col));
     }
   }
   distance = sqrt(distance);
   cout<<"distance="<<distance<<endl;
   if (distance <1.5){
     testError = 0;
   }
   else{
     testError = 1;
   }
   cout<<"testError="<<testError<<endl;
 }
 
 return testError;
 
}

