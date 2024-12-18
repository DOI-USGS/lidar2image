#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <iostream>
#include <cmath>

#include "ImageProcessing.h"

using namespace cv;
using namespace std;
using namespace at;

int main(int argc, char* argv[])
{
    // Open the file.
    Mat img = imread("test.jpg");
    if (img.empty()) 
      {
	cout << "Error: Couldn't open the image file." << endl;
	return 1;
      }
    
    //turn image into grayscale for easy testing
    Mat gray;
    cvtColor(img, gray, CV_BGR2GRAY);
    
    
    unsigned char* output_data1 = new unsigned char[img.rows * img.cols];
    float* output_data2 = new float[img.rows * img.cols];
    
    //try unsigned chars
    GaussianBlur<unsigned char>((unsigned char*)gray.data, gray.cols, gray.rows, output_data1, 7, 3);
    Mat blur_img1(gray.rows, gray.cols, CV_8UC1, output_data1);

    //try floats
    Mat float_img; 
    gray.convertTo(float_img, CV_32FC1);
    
    GaussianBlur<float>((float*)float_img.data, float_img.cols, float_img.rows, output_data2, 7, 3);
    Mat blur_img2(float_img.rows, float_img.cols, CV_32FC1, output_data2);
    
    //compare to opencv implementation
    Mat gray_blur_cv(gray.rows, gray.cols, CV_32FC1); 
    GaussianBlur(gray, gray_blur_cv, Size(7,7), 3, 3);
    gray_blur_cv.convertTo(gray_blur_cv, CV_32FC1);
    
    Mat diff(gray.rows, gray.cols, CV_32FC1);
    subtract(gray_blur_cv, blur_img2, diff);
    diff = abs(diff);
    
    cv::Scalar mean, mean_base, std_temp;
    cv::meanStdDev(diff, mean, std_temp);
    cv::meanStdDev(gray_blur_cv, mean_base, std_temp);
    
    cout << "mean difference is: " << mean/mean_base * 100 << "%" << endl;
    
    
    // Display the image.
    namedWindow("original", CV_GUI_EXPANDED);
    imshow("original", gray);
    
    namedWindow("blur1", CV_GUI_EXPANDED);
    imshow("blur1", blur_img1);
    
    namedWindow("blur2", CV_GUI_EXPANDED);
    imshow("blur2", blur_img2 / 255.0);

	namedWindow("opencv diff", CV_GUI_EXPANDED);
    imshow("opencv diff", diff/255);

	

    
    cvWaitKey(0);
    
    return 0;
}
