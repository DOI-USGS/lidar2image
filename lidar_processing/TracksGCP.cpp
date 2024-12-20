/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <iomanip>
#include <math.h>

#include "TracksLOLA.h"
#include "TracksFeatures.h"
#include "TracksGCP.h"
#include "../common/StringUtils.h"
//undefine thse macros in PCL to avoid conflicts with ISIS.
#undef RAD2DEG
#undef DEG2RAD

namespace at{

  vector<gcp> ComputeSalientLOLAFeatures(vector<vector< LOLAShot> > & LOLATracks)
  {
    int windowSize = 12;//16
    float salientFeatureThresh = 0.008;//0.015;
    
    // remove LOLA points that aren't discriminating
    vector<float> filter;
    filter =  MakeLidarFilter(windowSize);
    for (int ti = 0; ti < (int)LOLATracks.size(); ti++)
      ComputeSalientLOLAFeature(LOLATracks[ti], filter, salientFeatureThresh);
    
    vector<gcp> gcpArray;
    //this should be returned by ComputeSalientLOLAFeature
    for (unsigned int t=0; t<LOLATracks.size(); t++){
      for (unsigned int s=0; s<LOLATracks[t].size(); s++){
	if (LOLATracks[t][s].featurePtLOLA==1){
	  gcp this_gcp;
	  this_gcp.lon = LOLATracks[t][s].LOLAPts[2].x;
	  this_gcp.lat = LOLATracks[t][s].LOLAPts[2].y; 
	  this_gcp.rad = LOLATracks[t][s].LOLAPts[2].z*1000;
	  this_gcp.sigma_lon = 1.0;
	  this_gcp.sigma_lat = 1.0;
	  this_gcp.sigma_rad = 1.0;
	  this_gcp.year = LOLATracks[t][s].LOLAPts[2].year;
          this_gcp.month = LOLATracks[t][s].LOLAPts[2].month;
          this_gcp.day = LOLATracks[t][s].LOLAPts[2].day;
	  this_gcp.hour = LOLATracks[t][s].LOLAPts[2].hour;
	  this_gcp.min = LOLATracks[t][s].LOLAPts[2].min;
	  this_gcp.sec = LOLATracks[t][s].LOLAPts[2].sec;

	  this_gcp.trackIndex = t;
	  this_gcp.shotIndex = s;
	  gcpArray.push_back(this_gcp);
	}
      }
    }
    return gcpArray;
  }
  /*
  void SaveGCPoints( const vector<gcp>& gcpArray, const string& gcpFilename )
  {
    cout<<"num GCP arrays="<<gcpArray.size()<<endl;
    for( unsigned int i = 0; i < gcpArray.size(); ++i ){
      //check if this GCP is valid
      if( !gcpArray[i].filename.empty() ){
	stringstream filename;
	filename << gcpFilename << "_" << gcpArray[i].trackIndex 
		 << "_" << gcpArray[i].shotIndex << "_gcp.txt";
	
	ofstream file( filename.str().c_str() );
	if( !file ) {
	  cout<< "Can't open gcp output file " << filename.str()<<endl;
	}
	
	file << fixed << setprecision(6)
	     << gcpArray[i].lon << " "
	     << gcpArray[i].lat << " "
	     << gcpArray[i].rad << " "
	     << gcpArray[i].sigma_lon << " "
	     << gcpArray[i].sigma_lat << " "
	     << gcpArray[i].sigma_rad << " " << endl;
	
	for( unsigned int j = 0; j < gcpArray[i].x.size(); ++j ){
	  file << gcpArray[i].filename[j] << " " 
	       << fixed << setprecision(6) << gcpArray[i].x[j] << " " << gcpArray[i].y[j] << endl;
	}
	
	file.close();
      }
    }
  }
  */
  void SaveGCP( const vector<gcp>& gcpArray, const string& gcpFilename )
  {
    cout<<"num GCP arrays="<<gcpArray.size()<<endl;
    
    if ( !gcpFilename.empty() ){
      
      remove(gcpFilename.c_str());
     
      std::ofstream file;
      file.open(gcpFilename.c_str(), std::ios_base::app);
      
      if( !file ) {
	cout<< "Can't open gcp output file " << gcpFilename.c_str()<<endl;
      }
      
      for( unsigned int i = 0; i < gcpArray.size(); ++i ){
	
	file << fixed << setprecision(6)
	     << gcpArray[i].year<<"-"
             << ZeroPadNumber(gcpArray[i].month, 2)<<"-"
             << ZeroPadNumber(gcpArray[i].day, 2)<<"T" 
	     << ZeroPadNumber(gcpArray[i].hour, 2)<<":"
	     << ZeroPadNumber(gcpArray[i].min, 2)<<":"
	     << gcpArray[i].sec<< ","
	     << gcpArray[i].lon << ","
	     << gcpArray[i].lat << ","
	     << gcpArray[i].rad << ","
	     << gcpArray[i].sigma_lon << ","
	     << gcpArray[i].sigma_lat << ","
	     << gcpArray[i].sigma_rad << "," 
	     << GetFilenameNoPath(gcpArray[i].filename[0]) <<"," 
             << gcpArray[i].x[0] << "," 
             << gcpArray[i].y[0]<< endl;	  
      }
      
      file.close();
      
    }
  }
}
