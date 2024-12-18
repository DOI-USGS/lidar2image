// __BEGIN_LICENSE__
// Copyright (C) 2012 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#include <iostream>
#include <fstream>
#include "TracksLOLA.h"
#include "../common/FileListUtils.h"
#include "../common/StringUtils.h"

using namespace std;
using namespace at;


  


int main( int argc, char* argv[] )
{

  string lidarFilename;
  string resultsDirname="";
  
  int lidarFilenameFound = 0;
  int resultsDirnameFound = 0;
  int testFlag = 0;

  int testFailed = 0;
  
  for (int i = 1; i < argc; i++) { //We will iterate over argv[] to get the parameters stored inside.
   //Note that we're starting on 1 because we don't need to know the
   //path of the program, which is stored in argv[0]
   if (i + 1 != argc) { // Check that we haven't finished parsing already
     
     if (string(argv[i]).compare("-l")==0) {
       lidarFilename = argv[i + 1];
       lidarFilenameFound = 1;
     }
     else if (string(argv[i]).compare("-t")==0) {
       testFlag = 1;
     }    
     else if (string(argv[i]).compare("-r")==0) {
       resultsDirname = argv[i + 1];
       resultsDirnameFound = 1;
     }    
   }
 }

 int usageError = 0;
 
 if (lidarFilenameFound == 0){
   cout<<"Error! Lidar filename not found"<<endl;
   usageError = 1;
 }

 if (usageError == 1){
   std::cout << "Usage is "<< argv[0]<< " -l LOLAFilename -r resultsDirname"<<endl;
   std::cin.get();
   exit(0);
 }
 
 vector<string> outTrackFilenames;
 outTrackFilenames = trackSplitter(lidarFilename, resultsDirname);

 if (resultsDirname.compare("")!=0){
   string outTrackListFilename;
   outTrackListFilename = resultsDirname + "/"+GetFilenameNoPath(GetFilenameNoExt(lidarFilename))+"_list.txt";
   cout<<"outTrackListFilename="<<outTrackListFilename<<endl;
   ofstream outFile (outTrackListFilename.c_str());
   if (outFile.is_open()){
     for (int i=0; i<outTrackFilenames.size(); i++){
       outFile <<outTrackFilenames[i]<<endl;
     }
     outFile.close();
   }
   else cout << "Unable to open file "<<outTrackListFilename<<endl;
 }

 if (testFlag == 1){
   string testTrackListFilename = GetFilenameNoExt(lidarFilename)+"_test.txt";
   testFailed = test(outTrackFilenames, testTrackListFilename);
    
 }
 
 return testFailed;
}
