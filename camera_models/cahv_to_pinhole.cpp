#include <vector>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <iomanip>

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


#include "CAHV.h"
#include "pinhole.h"

using namespace std;
using namespace at;

int main( int argc, char *argv[] )
{
  if (argc<3){
    cout<<"Usage: cahv_to_pinhole cameraFilename resultsDirname"<<endl;
    exit(1);
  }

   string cameraFilename = string(argv[1]);
   string resultsDirname = string(argv[2]);
    
   cout << "Attempting to load file as a CAHV model" << endl;
   CAHV cahv(cameraFilename);
   
   cout << "Loaded CAHV model" << endl;
   
   // Construct a pinhole model from the cahv parameters to test the conversion
   double c[3], a[3], h[3], v[3];
   
   for(int i = 0; i < 3; i++){
     c[i] = cahv.m_c[i];
     a[i] = cahv.m_a[i];
     h[i] = cahv.m_h[i];
     v[i] = cahv.m_v[i];
   }
   
   int cahv_width = cahv.m_width;
   int cahv_height = cahv.m_height;
   
   // Create the pinhole
   cout << "Constructing a pinhole model" << endl;
   Pinhole converted(c, a, h, v, cahv_width, cahv_height);
   
   // Write the pinhole model to test
   string resultCamFilename = resultsDirname+"/pinhole.txt";
   converted.writeFile(resultCamFilename);
           
}







