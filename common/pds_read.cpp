#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "PdsToOpencv.h"

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

using namespace std;
using namespace cv;
using namespace at;

int main(int argc, char **argv )
{
  string filename;

  // Check the number of command line args:
  if (argc != 2) {
    cout << "Usage: \n" << "./pds_read <filename>" << endl; 
      exit(1);
  }
  else {
    filename = argv[1];
  }

  // Declarations for tile reader
  PDSTileReader tileReader; 
  long nCols, nRows, nBands, nBits;

  cv::Mat dstMat;

  IplImage *tile;

  // Initialize tile reader
  tileReader.initialize(filename);
  tileReader.getHeaderInfo(nCols, nRows, nBands, nBits);
  cout << "bands: " << nBands << " bits: " << nBits << endl;

  // Read tile to IplImage
  CvRect roi;
  roi.x = 0;
  roi.y = 0;
  roi.width = nCols/2;
  roi.height = nCols/2;


  int error = readPDSFileToMat(filename, dstMat);
  
  if (tileReader.readPDSTileToIplImage(tile, roi))
    exit(-1);
  
  cvSaveImage("test_tile.pgm", tile);
  cvReleaseImage(&tile);

  // Read file to IplImage
  IplImage *im;
  if (readPDSFileToIplImage(filename, im))
    exit(-1);

  cvSaveImage( "test_image.pgm", im);
  cvReleaseImage(&im);

  exit(0);
}



