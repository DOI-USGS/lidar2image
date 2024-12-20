/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/


//#include <omp.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "../common/StringUtils.h"
#include "PVLRead.h"

using namespace std;

namespace at
{
    
Eigen::Vector3f scanPVLForSunPosition(string pvlFilename)
{
  Eigen::Vector3f sunPosition;
//sunPosition.resize(3);

  ifstream fin (pvlFilename.c_str());
  string line;
  
  if (fin.is_open()){
    while (!fin.eof()){

      line.clear();
      std::getline(fin, line);
      FindAndReplace(line, ",","");
      FindAndReplace(line, "= (","");
 
      istringstream lineStream(line);
      string keyword("");
      lineStream >> keyword;
      
      if (keyword.compare(string("SunPosition"))==0){
	cout<<"line="<<line<<endl;
	lineStream>>sunPosition[0];
	lineStream>>sunPosition[1];

        line.clear();
	std::getline(fin, line);
	FindAndReplace(line, ") <km>","");
	cout<<"line="<<line<<endl;
	istringstream line2Stream(line);
	line2Stream>>sunPosition[2];
      }
    }
    cout<<"sunPosition="<<sunPosition[0]<<", "<<sunPosition[1]<<", "<<sunPosition[2]<<endl;
    fin.close();
  }
  else {
    cout<<"WARNING: could not open pvl file "+pvlFilename+"."<<endl;
  }
  return sunPosition;
}

Eigen::Vector3f scanPVLForCamPosition(string pvlFilename)
{
  Eigen::Vector3f spacecraftPosition;
  //spacecraftPosition.resize(3);

  ifstream fin (pvlFilename.c_str());
  string line;
  
  if (fin.is_open()){
    while (!fin.eof()){

      line.clear();
      std::getline(fin, line);
      FindAndReplace(line, ",","");
      FindAndReplace(line, "= (","");
 
      istringstream lineStream(line);
      string keyword("");
      lineStream >> keyword;
      
      if (keyword.compare(string("SpacecraftPosition"))==0){
	cout<<"line="<<line<<endl;
	lineStream>>spacecraftPosition[0];
	lineStream>>spacecraftPosition[1];

        line.clear();
	std::getline(fin, line);
	FindAndReplace(line, ") <km>","");
	cout<<"line="<<line<<endl;
	istringstream line2Stream(line);
	line2Stream>>spacecraftPosition[2];
      }
    }
    cout<<"spacecraftPosition="<<spacecraftPosition[0]<<", "<<spacecraftPosition[1]<<", "<<spacecraftPosition[2]<<endl;
    fin.close();
  }
  else {
    cout<<"WARNING: could not open pvl file "+pvlFilename+"."<<endl;
  }
  return spacecraftPosition;
}
                     
    
}
