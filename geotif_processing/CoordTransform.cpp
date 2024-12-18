#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include<iostream>

#include <gdal_priv.h>
#include <cpl_string.h>
#include <ogr_spatialref.h>

#include "../common/Tiling.h"
#include "CoordTransform.h"

using namespace std;

namespace at
{

    //from pixel to projected coordinates
    vector<double> pixelToProjected(vector<double> pixel, GDALDataset* data)
    {
        double geoTransform[6];
        data->GetGeoTransform(geoTransform);
        vector<double> projCoord;
        projCoord.resize(2);
        projCoord[0]=geoTransform[0] + pixel[0]*geoTransform[1] + pixel[1]*geoTransform[2];
        projCoord[1]=geoTransform[3] + pixel[0]*geoTransform[4] + pixel[1]*geoTransform[5];
        return (projCoord);
    }

    vector<double> projectedToPixel(vector<double> projCoord, GDALDataset* data)
    {
        double geoTransform[6];
        data->GetGeoTransform(geoTransform);
        vector<double> pixel;
        pixel.resize(2);
        pixel[1]=(geoTransform[1]*projCoord[1]-geoTransform[3]*geoTransform[1])/(geoTransform[1]*geoTransform[5]);
        pixel[0]=(projCoord[0]-geoTransform[0])/geoTransform[1];
        return pixel;
    }
    //from pixel to lon lat transformation
    vector<double> pixelToGeographic(vector<double> pixel, GDALDataset* data)
    {
      vector<double>geoCoord;
      double geoTransform[6];
      data->GetGeoTransform(geoTransform);
      cout.precision(12);
      
      OGRSpatialReference spatialRef;
      
      const char *wkt = data->GetProjectionRef();
      spatialRef.importFromWkt((char**)&wkt);
     
      if (spatialRef.IsProjected()){ 
	//cout<<"spatialRef is projected"<<endl;
	vector<double> projCoord;     
	projCoord.resize(3);
	projCoord[0]=geoTransform[0] + pixel[0]*geoTransform[1] + pixel[1]*geoTransform[2];
	projCoord[1]=geoTransform[3] + pixel[0]*geoTransform[4] + pixel[1]*geoTransform[5];
	//cout<<"geo="<<geoTransform[0]<<", "<<geoTransform[1]<<", "<<geoTransform[2]<<
	//", "<<geoTransform[3]<<", "<<geoTransform[4]<<", "<<geoTransform[5]<<endl;
	//cout<<"pixel="<<pixel[0]<<", "<<pixel[1]<<", frojCoord="<<projCoord[0]<<", "<<projCoord[1]<<endl;
	/*
	  vector<double> fProjCoord;     
	  fProjCoord.resize(3);
	  fProjCoord[0]=geoTransform[0] + pixel[0]*geoTransform[1] + pixel[1]*geoTransform[2];
	  fProjCoord[1]=geoTransform[3] + pixel[0]*geoTransform[4] + pixel[1]*geoTransform[5];
	  cout<<"pixel="<<pixel[0]<<", "<<pixel[1]<<", frojCoord="<<fProjCoord[0]<<", "<<fProjCoord[1]<<endl;
	  fProjCoord.resize(0);
	*/
	geoCoord = projectedToGeographic(projCoord, data);
	projCoord.resize(0);
	
      }
      else{
	cout<<"Warning: spatial ref is not projected"<<endl;
	geoCoord.resize(2);
	geoCoord[0]=geoTransform[0] + pixel[0]*geoTransform[1] + pixel[1]*geoTransform[2];
	geoCoord[1]=geoTransform[3] + pixel[0]*geoTransform[4] + pixel[1]*geoTransform[5];
      }  
      
      return geoCoord;
    }
  

    vector<double> geographicToPixel(vector<double> geoCoord, GDALDataset* data)
    {
        vector<double> pixel;
        pixel.resize(2);

        //cout<<"GEOGRAPHIC_TO_PIXEL"<<endl;
        OGRSpatialReference spatialRef;
        const char *wkt = data->GetProjectionRef();
        spatialRef.importFromWkt((char**)&wkt);
        //delete wkt;

        double geoTransform[6];
        data->GetGeoTransform(geoTransform);
        /*
          if (!spatialRef.IsProjected()){
          cout<<"Warning! spatial ref does not contain projection info. assume equirectagular projection"<<endl; 
          spatialRef.SetProjCS( "UTM 10 (WGS84) in northern hemisphere." );
          spatialRef.SetWellKnownGeogCS( "WGS84" );
          spatialRef.SetUTM( 10, TRUE );
          }
        */
        //geographic to projected
        vector<double> projCoord = geographicToProjected(geoCoord, data);
        //projected to pixel
        pixel[0]=(projCoord[0]-geoTransform[0])/geoTransform[1];
        pixel[1]=(projCoord[1]-geoTransform[3])/geoTransform[5];
  
        projCoord.resize(0);
  
        return pixel;
    }

    //projected (point) to geographic (lonlat) coordinates 
    vector<double> projectedToGeographic(vector<double> projectedCoord, GDALDataset* data)
    {

        vector<double> geographicCoord;
    
        OGRSpatialReference  *poLatLong;
        OGRCoordinateTransformation *poTransform;
    
        OGRSpatialReference spatialRef;
        const char *wkt = data->GetProjectionRef();
        spatialRef.importFromWkt((char**)&wkt);
        if (!spatialRef.IsProjected()){
            //cout<<"Warning! spatial ref does not contain projection info. assume equirectagular projection"<<endl; 
            spatialRef.SetProjCS( "UTM 10 (WGS84) in northern hemisphere." );
            spatialRef.SetWellKnownGeogCS( "WGS84" );
            spatialRef.SetUTM( 10, TRUE );
        }
    
        poLatLong = spatialRef.CloneGeogCS();
        poTransform = OGRCreateCoordinateTransformation( &spatialRef, poLatLong );
        if ( poTransform == NULL ){
            cerr << "ERROR [projectedToGeographic()]: OGRCoordinateTransformation could not be created" << endl;
        }
    
        int nPoints = 1;
        double *x=new double[nPoints];
        double *y=new double[nPoints];
        double *z=new double[nPoints];
        x[0] = projectedCoord[0];
        y[0] = projectedCoord[1];
        z[0] = projectedCoord[2];

        if ( !poTransform->Transform( nPoints, x, y, z ) ){
            cerr << "ERROR [projectedToGeographic()]: OGRCoordinateTransformation::Transform() failed!" << endl;
        }
  
        geographicCoord.resize(3);
        geographicCoord[0]=x[0];
        geographicCoord[1]=y[0];
        geographicCoord[2]=z[0];
        delete x; delete y; delete z;

        delete poTransform;
        delete poLatLong;
    
        return geographicCoord;
    }

    //projected (point) to geographic (lonlat) coordinates 
    vector<double> geographicToProjected(vector<double> geographicCoord, GDALDataset* data)
    {
        vector<double> projectedCoord;
    
        OGRSpatialReference  *poLatLong;
        OGRCoordinateTransformation *poTransform;
    
        OGRSpatialReference spatialRef;
        const char *wkt = data->GetProjectionRef();
        spatialRef.importFromWkt((char**)&wkt);

        if (!spatialRef.IsProjected()){
            cout<<"Warning! spatial ref does not contain projection info. assume equirectagular projection"<<endl; 
            spatialRef.SetProjCS( "UTM 10 (WGS84) in northern hemisphere." );
            spatialRef.SetWellKnownGeogCS( "WGS84" );
            spatialRef.SetUTM( 10, TRUE );
        }

        poLatLong = spatialRef.CloneGeogCS();
    
        poTransform = OGRCreateCoordinateTransformation( poLatLong, &spatialRef);
        if ( poTransform == NULL ){
            cerr << "ERROR [geographicToProjected()]: OGRCoordinateTransformation could not be created" << endl;
        }

        int nPoints = 1;
        double *x=new double[nPoints];
        double *y=new double[nPoints];
        double *z=new double[nPoints];
        x[0] = geographicCoord[0];
        y[0] = geographicCoord[1];
        z[0] = geographicCoord[2];

        if ( !poTransform->Transform( nPoints, x, y, z ) ){
            cerr << "ERROR [geographicToProjected()]: OGRCoordinateTransformation::Transform() failed!" << endl;
        }
        projectedCoord.resize(3);
        projectedCoord[0]=x[0];
        projectedCoord[1]=y[0];
        projectedCoord[2]=z[0];

        delete x; delete y; delete z;
        delete poLatLong;
        delete poTransform;

        return projectedCoord;
    }

    std::vector<double> geographicToCartesian(std::vector<double> geoCoord, GDALDataset* data)
    {
      std::vector<double> cartesianCoord;

      const char *wkt = data->GetProjectionRef();
      OGRSpatialReferenceH hSRS((char**)&wkt);
  
      //Get spheroid semi major axis.
      double semiMajor = OSRGetSemiMajor(hSRS, NULL);
      //Get spheroid semi minor axis. 			
      double semiMinor = OSRGetSemiMinor(hSRS, NULL);
      cout<<"semiMinor = "<<semiMinor<<", semiMajor="<<semiMajor<<endl;
 			
      return cartesianCoord;
    }

    Eigen::Vector3f geographicToCartesian(Eigen::Vector3f lonLatRad)
    {

      Eigen::Vector3f cartesianPt;
      //float radius = 1737400.0;
      //compute global coordinates for the height map.
      
      double lonRad   = lonLatRad[0]; 
      double latRad   = lonLatRad[1];
      double altitude = lonLatRad[2];
      
   
      cartesianPt(0) = (altitude) * cos(latRad) * cos(lonRad);
      cartesianPt(1) = (altitude) * cos(latRad) * sin(lonRad);
      cartesianPt(2) = (altitude) * sin(latRad);
      
      return cartesianPt;
    }
  
    Eigen::Vector3f geographicToCartesian(Eigen::Vector3f lonLatRad, double radius)
    {
      double lonRad  = lonLatRad[0]*(M_PI/180.0);
      double latRad  = lonLatRad[1]*(M_PI/180.0);
      double altitude= lonLatRad[2];
      
      Eigen::Vector3f cartesianPt;
      cartesianPt(0) = (radius + altitude) * cos(latRad) * cos(lonRad);
      cartesianPt(1) = (radius + altitude) * cos(latRad) * sin(lonRad);
      cartesianPt(2) = (radius + altitude) * sin(latRad);
      
      return cartesianPt;
    }
  
  

}
