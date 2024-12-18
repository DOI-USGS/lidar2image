#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Tiling.h"
#include "FileListUtils.h"
using namespace std;
using namespace at;


void printUsage()
{
  cout << endl;
  cout << "****************** USAGE *****************" << endl;
  cout << "image_tiling <configFile> <inputFile>              " << endl;
  cout << "******************************************" << endl;
  cout << endl;
}

int main(int argc, char** argv)
{
  vector<string> inputFiles;
  Tiling tiler;
  string inputFilename;
  string outputFilename;
  string configFilename;
  int frameIndex;
  cv::Mat image;
  
  //Read Command Line Arguments
  if(argc != 3){
    printUsage();
    return (-1);
  }
  
  configFilename = string(argv[1]);
  inputFilename = string(argv[2]);
  
  //Read Tile Configuration File
  tiler = Tiling(configFilename);
  
  ReadFileList(inputFilename, inputFiles);

  for(frameIndex=0;frameIndex<inputFiles.size();frameIndex++){
    
    cout<<"inputFilename="<<inputFiles[frameIndex]<<endl;
    image = cvLoadImage(inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
    
    tiler.process(image.rows, image.cols);
    
    for(int i=0; i<tiler.m_tiles.size(); i++){
      //create a CvRect from the bbox structure
      cv::Rect thisRect;
      thisRect.x = tiler.m_tiles[i].xl;
      thisRect.y = tiler.m_tiles[i].yt;
      thisRect.width = tiler.m_tiles[i].width;
      thisRect.height = tiler.m_tiles[i].height;
      
      cout << "BB " << i << ": [" << thisRect.x << ", " << thisRect.y 
	   <<  ", " << thisRect.width << ", " << thisRect.height << "]" << endl; 
      
      cv::Mat tileImage = image(thisRect);
      imshow("Tile", tileImage);
      
      cvWaitKey(1000);
    }
  }
  
  return 0;
}

