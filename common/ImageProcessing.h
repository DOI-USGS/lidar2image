#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_

#include <vector>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;

namespace at
{
      
  bool write_pgm_byte_image(uint8_t *pBuffer, const char* filename,
			    const long nCols, const long nRows);
  bool write_ppm_byte_image(uint8_t *rgbBuffer, const char* filename,
			    const long nCols, const long nRows,
			    bool bandSequential = true);
  unsigned char* read_pgm(std::string filename, long *nCols, long *nRows);
  bool writeFloatToBytePGM(float *data, string filename, const long imageWidth,
			   const long imageHeight, float gain, float bias);

  //simple hole filling method for linear 1d holes uses noDataValue or NaN values to determine holes
  void holeFilling(float *dem, int width, int height, float noDataValue);
  //computes the bilinear interpolation function
  float GetInterpValue(float x, float y, float *image, int width, int height);
  
  // function GaussianBlur_noCV - symmetric gaussian blur of an image using two pass approach, does not use OpenCV
  // It operates on a single channel i.e. NxMx1 gray image, so you will need to call 3x for a color image
  // Data is assumed to be row ordered: i.e. im(x,y) = im[y*row_size + x]
  // A "replicates" border functionality is implemented (takes the nearest value in convolution) 
  // 
  // x_width: i.e. image.cols
  // y_height: i.e. image.rows
  // output_img: new unsigned char [image.cols * image.rows] - you manage this memory
  // mask_size: the width/height of the square mask
  // sigma: the standard deviation of the gaussian kernel
  //
  // return false if failed, true if success
  template <typename StorageType> bool
    GaussianBlur(StorageType* input_img, int x_width, int y_height, StorageType* output_img, int mask_size, double sigma)
  {
    // mask dimension must be odd
    if (!(mask_size % 2))
      return false;
    
    if (!output_img || !input_img)
      return false;
    
    int pivot = mask_size/2;
    
    //generate lookup table for mask
    double* gauss_mask = new double[mask_size];
    //double scalar = 1.0 / (sigma * sqrt(2 * 3.14159)); don't need, this will cancel out anyways
    double normalizer = 0.0;
    for (int i = 0; i < mask_size; i++)
    {
	double exponent = exp(-(double)((i-pivot)*(i-pivot))  / (2.0 * sigma * sigma));
	gauss_mask[i] = exponent;
	normalizer += exponent; 
    }
    
    //normalize to unit energy
    for (int i = 0; i < mask_size; i++)
      {
	gauss_mask[i] /= normalizer;
	cout << (double)gauss_mask[i] << " ";
      }
    cout << endl << endl;
    
    //create temporary image
    StorageType* temp_buffer = new StorageType[x_width * y_height];
    
    //horizontal pass
    for (int j = 0; j < y_height; j++)
      {
	for (int i = 0; i < x_width; i++)
	  {
	    double conv_sum = 0;
	    
	    //calculate the local convolution
	    for (int nhood = -pivot; nhood <= pivot; nhood++)
	      {
		//check for boundary conditions
		//in the case of the boundary, the action is to 
		//"replicate" the nearest value
		int loc = i+nhood; 
		if (loc < 0) 
		  loc = 0;
		else if (loc >= x_width)
		  loc = x_width -1;
		
		//add the contribution of this neighbor 
		conv_sum += gauss_mask[nhood + pivot] * (double)input_img[j*x_width + loc]; 
	      }
	    
	    temp_buffer[j*x_width + i]= (StorageType)conv_sum;
	    //temp_buffer[j*x_width + i] = input_img[j*x_width + i];	//pass thru
	  }
      }
    
    //vertical pass
    for (int i = 0; i < x_width; i++)
      {
	for (int j = 0; j < y_height; j++)
	  {
	    double conv_sum = 0;
	    
	    //calculate the local convolution
	    for (int nhood = -pivot; nhood <= pivot; nhood++)
	      {
		//check for boundary conditions
		//in the case of the boundary, the action is to 
		//"replicate" the nearest value
		int loc = j+nhood; 
		if (loc < 0) 
		  loc = 0;
		else if (loc >= y_height)
		  loc = y_height -1;
		
		//add the contribution of this neighbor 
		conv_sum += gauss_mask[nhood + pivot] * (double)temp_buffer[loc*x_width + i]; 
	      }
	    
	    output_img[j*x_width + i] = (StorageType)conv_sum;
	    //output_img[j*x_width + i] = temp_buffer[j*x_width + i];	//pass thru
	  }
      }
    
    //delete memory
    delete [] temp_buffer;
    delete [] gauss_mask;
    
    return true;
  }

  //sun angle adjustment or in other words albedo reconstruction, uses lambertian model, pretty cool
  template <typename ValueType>
    void
    sun_angle_adjustment(ValueType *valueBuffer, const uint32_t numValues,
                         const float sunElevation)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      
      // Get IOF image lighting scale factor
      float iofScale = fabs(cos((90.0 + sunElevation)*M_PI/180.0));
      // We don't let iofScale get smaller than 0.1 (correspond to about 85 deg
      // elevation, because we're dividing current pixel value with it and
      // don't want to grow too huge.
      if (iofScale < 0.1)
	iofScale = 0.1;
      iofScale = 1.0/iofScale;
      
      for (uint32_t i = 0; i < numValues; i++){
	double value = iofScale * valueBuffer[i];
	valueBuffer[i] = (value > kMaxTypeValue) ? kMaxTypeValue : ValueType(value);
      }
    }
  
  
  //all values are multiplied by a constant gain and summed to a constant bias
  template <typename ValueType>
    void
    radiance_calibration_adjustment(ValueType *valueBuffer, const uint32_t numValues,
                                    const double radianceScale,
                                    const double radianceOffset)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      
      for (uint32_t i = 0; i < numValues; i++){
	double value = radianceScale * valueBuffer[i] + radianceOffset;
	valueBuffer[i] = (value > kMaxTypeValue) ? kMaxTypeValue : ValueType(value);
      }
    }
  
  
  //all values are multiplied by a fixed contant
  template <typename ValueType>
    void
    exposure_adjustment(ValueType *valueBuffer, const uint32_t numValues,
                        const double exposureRatio)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      
      for (uint32_t i = 0; i < numValues; i++){
	double value = exposureRatio * valueBuffer[i];
	valueBuffer[i] = (value > kMaxTypeValue) ? kMaxTypeValue : ValueType(value);
      }
    }
  
  template <typename ValueType>
    uint32_t *
    generate_histogram(ValueType valueBuffer[], const uint32_t numValues,
                       ValueType &minVal, ValueType &maxVal)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      
      // compute min/max and image histogram
      minVal = kMaxTypeValue;
      maxVal = 0;
      uint32_t numPixelValues = kMaxTypeValue + 1;
      uint32_t *histogram = new uint32_t[numPixelValues];
      
      for (uint32_t i = 0; i < numValues; i++){
	const uint16_t pixval = valueBuffer[i];
	
	histogram[pixval] += 1;
	
	if (pixval < minVal)
	  minVal = pixval;
	if (pixval > maxVal)
	  maxVal = pixval;
      }
      
      cout << "generate_histogram(): input image (min, max) = (" << double(minVal)
	   << ", " << double(maxVal) << ")" << endl;
      
      // Now, we count up and throw out the outliers at the top end of the
      // image histogram...
#if 0
      int totalPixels = 0;
      for (int i = maxValue; i >= 0; i--){
	totalPixels += histogram[i];
	if (totalPixels > (HISTOGRAM_OUTLIER_PERCENTAGE * (float)numValues))
	  {
	    maxIndex = i-1;
	    break;
	  }
      }
      printf("Adjusted max: %1d, Min: %1d\n", maxVal, minVal);
#endif
      
      return histogram;
    }
  
  //all values are scaled such that the max value is constant. 
  template <typename ValueType>
    void
    auto_gain_adjustment(ValueType *valueBuffer, const uint32_t numValues)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      ValueType min, max;
      
      generate_histogram(valueBuffer, numValues, min, max);
      double autoGain = double(kMaxTypeValue) / double(max);
      
      cout << "auto_gain_adjustment() max type value = " << double(kMaxTypeValue) << endl;
      cout << "auto_gain_adjustment() calculated gain = " << autoGain << endl;
      
      for (uint32_t i = 0; i < numValues; i++){
	double value = double(valueBuffer[i]) * autoGain;
	valueBuffer[i] = ((value > double(kMaxTypeValue)) ?
			  kMaxTypeValue : ValueType(round(value)));
      }
    }
  
  template <typename ValueType>
    void
    user_gain_adjustment(ValueType *valueBuffer, const uint32_t numValues,
                         double userScale)
    {
      const static ValueType kMaxTypeValue = numeric_limits<ValueType>::max();
      
      for (uint32_t i = 0; i < numValues; i++){
	double value = double(valueBuffer[i]) * userScale;
	valueBuffer[i] = (value > kMaxTypeValue) ? kMaxTypeValue : ValueType(value);
      }
    }
    
  //****************Color Norm functions**************** 
  template <typename ChannelType> 
    static inline ChannelType colorL1Norm(std::vector<ChannelType> & color)
  {
    if (numeric_limits<ChannelType>::is_integer) {
      if (numeric_limits<ChannelType>::is_signed){
	return (abs(color[0]) + abs(color[1]) + abs(color[2]));
      }
      else {
	//since this is a distance metric, if an unsigned type overflows, 
	//we do not want it to wrap around, e.g. 255+255+255 = 253 for uint8.
	//This code will saturate at the max value for that data type.
	double out = color[0] + color[1] + color[2];
	out = (out > numeric_limits<ChannelType>::max() ) ? numeric_limits<ChannelType>::max() : out;
	
	return out;
      }
    }
    else {
      return (fabs(color[0]) + fabs(color[1]) + fabs(color[2]));
    }
  }
  
  template <typename ChannelType> 
    static inline double colorL2Norm(std::vector<ChannelType> & color)
  {
    return sqrt(color[0]*color[0] +
		color[1]*color[1] +
		color[2]*color[2]);
  }

  template <typename ChannelType> 
    static inline double colorL2NormSqr(std::vector<ChannelType> & color)
  {
    return (color[0]*color[0] +
	    color[1]*color[1] +
	    color[2]*color[2]);
  }

  template <typename ChannelType>
    static inline double colorLNNorm(std::vector<ChannelType> & color, double N)
  {
    return pow(pow(color[0], N) +
	       pow(color[1], N) +
	       pow(color[2], N),
	       1.0/N);
  }

  template <typename ChannelType> 
    static inline double colorLInfNorm(std::vector<ChannelType> & color)
  {
    double maxTemp = ((color[0] > color[1]) ? color[0] :
		      color[1]);
    return  (maxTemp > color[2]) ? maxTemp : color[2];
  }

  template <typename ChannelType>
    static inline double colorAngle(std::vector<ChannelType> & colorA,
	       std::vector<ChannelType> & colorB)
  {
    double dotProduct = 0;
    
    for (int i = 0; i < 3; i++)
      dotProduct += (double)colorA[i] * (double)colorB[i];
    
    return acos(dotProduct/(colorL2Norm(colorA)*colorL2Norm(colorB)));
  }
}

#endif	// PDS_IMAGE_PROCESSING_H_
