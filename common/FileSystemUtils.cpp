#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>

//
#include <vector>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <unistd.h>
#include <stdlib.h>

#include "FileSystemUtils.h"
#include "StringUtils.h"

using namespace std;
 
namespace at
{
  
  bool is_file_exist(const char *fileName)
  {
    std::ifstream infile(fileName);
    return infile.good();
  }


  //returns in variable files the filenames in directory dir
  //Names are returned in an arbitrary order, sorting is recommended
  int getFileNames (string dir, vector<string> &files)
  {
    string fileName;
    DIR *dp;
    struct dirent *dirp;
    string name;
    
    files.clear();
    
    //open directory
    if((dp  = opendir(dir.c_str())) == NULL) {
      cout << "Error opening " << dir << endl;
      return -1;
    }
    
    //read file names
    while ((dirp = readdir(dp)) != NULL) {
      name=dirp->d_name;
      fileName=dir+'/'+name;
      std::string ext = GetFilenameExt(fileName);
      if (ext.compare(string("IMG"))==0){
	files.push_back(fileName);
      }
    }
    
    closedir(dp);
    return 0;
  }
  
  //returns in variable files the filenames in directory dir
  //Names are returned in an arbitrary order, sorting is recommended
  int GetImageFileNames (string dir, vector<string> &files)
  {
    string fileName;
    DIR *dp;
    struct dirent *dirp;
    string name;
    
    files.clear();
    
    //open directory
    if((dp  = opendir(dir.c_str())) == NULL) {
      cout << "Error opening " << dir << endl;
      return -1;
    }
    
    //read file names
    while ((dirp = readdir(dp)) != NULL) {
      name=dirp->d_name;
      fileName=dir+'/'+name;
      std::string ext = GetFilenameExt(fileName);
      if ((ext.compare(string("IMG"))==0) || (ext.compare(string("JPG"))==0) || (ext.compare(string("jpg"))==0) ||
	  (ext.compare(string("bmp"))==0) || (ext.compare(string("pgm"))==0) || (ext.compare(string("png"))==0) ||
	  (ext.compare(string("tiff"))==0)|| (ext.compare(string("tif"))==0)){
	files.push_back(fileName);
      }
    
    }
    
    closedir(dp);
    return 0;
  }
  
  // Returns a list of all files, regardless of extension
  // includeDirs will include directory names from results if set to true
  // Names are returned in an arbitrary order, sorting is recommended
  int getFileNames (string dir, vector<string> &files, bool includeDirs)
  {
    string fileName;
    DIR *dp;
    struct dirent *dirp;
    string name;

    struct stat sb;
    
    files.clear();
    
    //open directory
    if((dp  = opendir(dir.c_str())) == NULL) {
      cout << "Error opening " << dir << endl;
      return -1;
    }
    
    //read file names
    while ((dirp = readdir(dp)) != NULL) {
      name=dirp->d_name;
      fileName=dir+'/'+name;

      // Check if this is a file and not a directory
      // Exclude directory names if includeDirs is set to false
      if (includeDirs || !(stat(fileName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))){
	files.push_back(fileName);
      }
    }
    
    closedir(dp);
    return 0;
  }

 
}


