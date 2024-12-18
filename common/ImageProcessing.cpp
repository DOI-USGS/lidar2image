#include <iostream>
#include <limits>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <iostream> // cout, cerr
#include <fstream> // ifstream
#include <sstream> // stringstream

#include "ImageProcessing.h"
using namespace std;

namespace at
{

    /////////////////////////////////////////////////////////////////////////////
    // Writes PGM byte image

    bool write_pgm_byte_image(uint8_t *pBuffer, const char* filename,
                              const long nCols, const long nRows)
    {
        FILE* file = fopen(filename, "w");
        if (!file) return(false);
        
        // write magic number (for binary data)
        fprintf(file, "P5\n");
        // write width, height
        fprintf(file, "%ld %ld\n", nCols, nRows);
        // write max grayscale value
        fprintf(file,"255\n"); // max color val

        // write the image in binary format to file
        const long size = nCols * nRows * sizeof(uint8_t);
        if (!fwrite(pBuffer, 1, size, file))
            return false;

        fclose(file);

        return true;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Writes a PPM byte image

    bool write_ppm_byte_image(uint8_t *rgbBuffer, const char* filename,
                              const long nCols, const long nRows,
                              bool bandSequential)
    {
        FILE* outFP = fopen(filename, "wb");
        if (!outFP)
            return(false);
        
        // write magic number (for binary data)
        fprintf(outFP, "P6\n");
        // write width, height
        fprintf(outFP, "%ld %ld\n", nCols, nRows);
        // write max grayscale value
        fprintf(outFP,"255\n"); // max color val

        const long numPixels = nRows * nCols;
        // write the image in binary format to file
        if (bandSequential)		// Image bytes are RGBRGBRGB...
        {
            const long bufferSize = numPixels * 3 * sizeof(uint8_t);
            if (!fwrite(rgbBuffer, 1, bufferSize, outFP)) 
                return false;
        }
        else		     // Image has separate red, green, and blue planes
        {
            const long planeSize = numPixels * sizeof(uint8_t);
            const long rgbRowBufferSize = 3 * nCols * sizeof(uint8_t);
            uint8_t *redBuffer = rgbBuffer;
            uint8_t *greenBuffer = redBuffer + planeSize;
            uint8_t *blueBuffer = greenBuffer + planeSize;
            uint8_t rgbRowBuffer[rgbRowBufferSize];
            long planeIndex = 0;
            for (long row = 0; row < nRows; ++row)
            {
                // For speed, we write a row at a time (i.e., rather than one
                // pixel at at time)
                long rowIndex = 0;
                for (long col = 0; col < nCols; ++col)
                {
                    rgbRowBuffer[rowIndex] = redBuffer[planeIndex];
                    ++rowIndex;
                    rgbRowBuffer[rowIndex] = greenBuffer[planeIndex];
                    ++rowIndex;
                    rgbRowBuffer[rowIndex] = blueBuffer[planeIndex];
                    ++rowIndex;
                    ++planeIndex;
                }
                fwrite(rgbRowBuffer, 1, rgbRowBufferSize, outFP);
            }
        }

        fclose(outFP);

        return true;
    }

    bool writeFloatToBytePGM(float *data, string filename, const long imageWidth,
			     const long imageHeight, float gain, float bias)
    {
      
      unsigned char*scaleImageData = new unsigned char[imageWidth*imageHeight];
      for (int i = 0; i < imageWidth*imageHeight; i++){
	scaleImageData[i] = (unsigned char)(round(gain*data[i]+bias));
      }
      
      bool flag = write_pgm_byte_image(scaleImageData,
				       filename.c_str(),
				       imageWidth,
				       imageHeight);
      delete scaleImageData;
      
      return flag;
      
    }
  
    unsigned char* read_pgm(string filename, long *numCols, long *numRows)
    {

      FILE* file = fopen(filename.c_str(), "r");
      if (!file) return(nullptr);
      
      char *format; format=new char[2];
      fscanf(file, "%s\n", format);
      cout<<"format="<<format<<endl;
      long numrows, numcols;
      fscanf(file, "%ld %ld\n", &numcols, &numrows);
      cout<<"numCols="<<numcols<<", numRows="<<numrows<<endl;
      int maxValFile;
      fscanf(file,"%d\n", &maxValFile); // max color val
      cout<<"maxVal="<<maxValFile<<endl;
      unsigned char *image = new unsigned char[numrows*numcols];
      size_t result = fread (image, 1, numcols*numrows, file);
   
      fclose(file);

    
      *numRows = numrows;
      *numCols = numcols;
      return image;
    }
  //simple hole filling method for linear 1d holes uses noDataValue or NaN values to determine holes
  void holeFilling(float *dem, int width, int height, float noDataValue)
  {

    int rowStartIndex, index, holeStart, holeEnd;
 
    for (int row = 0; row < height; row++){
      rowStartIndex = row*width;
      if ((dem[rowStartIndex]==noDataValue)||(isnan(dem[rowStartIndex]))){
	holeStart = rowStartIndex;
      }
      for (int col = 1; col < width; col++){
	index = rowStartIndex + col;
	
	if ((dem[index]==noDataValue || isnan(dem[index])) && (dem[index-1]!=noDataValue && !isnan(dem[index-1]))){
	  holeStart=col;//first invalid value
	}
	
	if ((dem[index]!=noDataValue && !isnan(dem[index])) && (dem[index-1]==noDataValue || isnan(dem[index-1]))){
	  holeEnd=col-1;//last invalid value
	  //std::cout<<" rowstart = "<<rowStartIndex<<", start = "<<holeStart<<", end="<<holeEnd<<", diff = "<<holeEnd-holeStart<<std::endl;
	  if (holeStart!=rowStartIndex){
	    //fill the hole from holeStart to holeEnd
	    float alpha = (dem[rowStartIndex+holeEnd+1]-dem[rowStartIndex+holeStart-1])/(holeEnd-holeStart+1);
	    for (int k = holeStart; k <= holeEnd; k++){
	      dem[rowStartIndex+k]=(k-holeStart+1)*alpha+dem[rowStartIndex+holeStart-1];
	    }
	  }
	}
	
      }
    }
  }

   float GetInterpValue(float x, float y, float *image, int width, int height)
  {

    float interp;
    int lx, rx, ty, by;
    lx = floor(x);
    ty = floor(y); 
    rx = lx+1;
    by = ty+1;
  
    //conditions to guarantee that the neighborng pixels 
    //are inside the image boundaries
    if (lx < 0){
        lx = 0;
    }
    if (lx > width-1){
        lx = width-1;
    }
    
    if (rx < 0){
        rx = 0;
    }
    if (rx > width-1){
        rx = width-1;
    }

    if (ty < 0){
        ty = 0;
    }
    if (ty > height-1){
        ty = height-1;
    }
    
    if (by < 0){
        by = 0;
    }
    if (by > height-1){
        by = height-1;
    }
    
    int leftTop = ty*width+lx;
    int rightTop = leftTop+1;
    if (rightTop > width*height-1){
        rightTop = width*height-1;
    }
    int leftBottom = leftTop+width;
    if (leftBottom > width*height-1){
        leftBottom = width*height-1;
    }
    int rightBottom = leftBottom+1;
    if (rightBottom > width*height-1){
        rightBottom = width*height-1;
    } 
    float tl = image[ty*width+lx];//image.at<float>(ty, lx);
    float bl = image[by*width+lx];//image.at<float>(by, lx);
    float tr = image[ty*width+rx];//image.at<float>(ty, rx);
    float br = image[by*width+rx];//image.at<float>(by, rx);
    
    interp = tl*(rx-x)*(by-y)+
      tr*(x-lx)*(by-y)+
      bl*(rx-x)*(y-ty)+
      br*(x-lx)*(y-ty);
    
    interp = interp/((rx-lx)*(by-ty));
 
    return interp;
  };

}
