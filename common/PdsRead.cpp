// $Id$
// 
// pds2pgm: convert an unsigned short PDS image into a pgm image.
// 
// This program also reads from the input PDS image's header the CAHV
// camera parameters, translates them into a pan/tilt pair of values,
// and writes those into a comment in the output pgm image's header,
// that Viz can understand.
//
// Author: Doron Tal
// Date Created: Dec. 11, 2003

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

extern "C" { 
#include "oal.h"
#include "lablib3.h"
}

#include "PdsRead.h"

using namespace std;

namespace at
{

    long get_pds_label_value(ODLTREE odltree, const char* strLabelName) 
    {
        long result;
        OaKwdValuetoLong((char*)strLabelName, odltree, &result);
        if (result == 0) {
            cerr << "get_pds_label_value(): Error reading label "
                 << strLabelName << endl;
            return -1;
        }
        return result;
    }

    /////////////////////////////////////////////////////////////////////////////
    int getRadianceDataForImage(ODLTREE odltree, double *radianceScaleFactor,
                                double *radianceOffset)
    {
        KEYWORD *groupKeyword, *parameterKeyword;
        groupKeyword = OdlFindKwd(odltree, (char*) "GROUP", (char*) "DERIVED_IMAGE_PARMS", 1, 
                                  ODL_THIS_OBJECT);
        if (groupKeyword == NULL) {
            cerr << "getRadianceDataForImage(): Could not find DERIVED_IMAGE_PARMS!"
                 << endl;
            return -1;
        }
  
        parameterKeyword = OdlNextKwd( groupKeyword, (char*) "RADIANCE_SCALING_FACTOR", (char*) "*",
                                       1, ODL_THIS_OBJECT);
        if (parameterKeyword == NULL) {
            cerr << "getRadianceDataForImage(): Could not find RADIANCE_SCALING_FACTOR"
                 << endl;
            return -1;
        }

        char *conversionPos = 0;
        char *keywordValue = OdlGetKwdValue(parameterKeyword);

        if (strcmp(keywordValue, "\"NULL\"") == 0)
            return -1;

        *radianceScaleFactor = strtod(keywordValue, &conversionPos);
        if (conversionPos == keywordValue)
        {
            cerr << "getRadianceDataForImage(): invalid RADIANCE_SCALING_FACTOR value"
                 << endl;
            return -1;
        }

        parameterKeyword = OdlNextKwd(groupKeyword, (char*) "RADIANCE_OFFSET", (char*) "*",
                                      1, ODL_THIS_OBJECT);
        if (parameterKeyword == NULL) {
            cerr << "getRadianceDataForImage(): Could not find RADIANCE_OFFSET!"
                 << endl;
            return -1;
        }

        keywordValue = OdlGetKwdValue(parameterKeyword);
        if (strcmp(keywordValue, "\"NULL\"") == 0)
            return -1;

        *radianceOffset = strtod(keywordValue, &conversionPos);
        if (conversionPos == keywordValue)
        {
            cerr << "getRadianceDataForImage(): invalid RADIANCE_OFFSET value"
                 << endl;
            return -1;
        }

        return 0;
    }
    int getRadianceDataForImage(string strInFilename, double *radianceScaleFactor,
                                double *radianceOffset)
    {
        // parse the pds file to get the top level node odltree
        ODLTREE odltree = OaParseLabelFile((char *)strInFilename.c_str(), 
                                           (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);
  
        if (odltree == NULL) { 
            OaReportError((char*) "Error in OaParseLabelFile()");
            return -1; 
        }

        getRadianceDataForImage(odltree, radianceScaleFactor,
                                radianceOffset);

        return 0;
    }
    /////////////////////////////////////////////////////////////////////////////

    int scan_file_for_string(FILE *fp, const char *str)
    {
        const char *s;
        int c;
    
        s = str;
        for (;;) {
            c = getc(fp);
            if (c == EOF)
                return -1;
            if (c == *s) {
                s++;
                if (*s == '\0')
                    return 0;
            }
            else
                s = str;
        }
    }

    /////////////////////////////////////////////////////////////////////////////

    int scan_file_for_cahv_params(const char* strInFilename, 
                                  float c[3], float a[3], float h[3], float v[3],
                                  float o[3], float r[3], bool &gotOR)
    {
        gotOR = true;
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // first scan for the camera parameters group
        if (scan_file_for_string(hFile, CAM_PARAMS_MARKER)) {
            cerr << " error: could not find camera parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
    
        // scan and grab center axis
        if (scan_file_for_string(hFile, C_MARKER)) {
            cerr << " error: could not find camera parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        } 
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &c[0], &c[1], &c[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

#ifdef _DEBUG
            cerr << "C: " << c[0] << " " << c[1] << " " << c[2] << endl;
#endif
        }

        if (scan_file_for_string(hFile, A_MARKER)) {
            cerr << " error: could not find camera parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &a[0], &a[1], &a[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "A: " << a[0] << " " << a[1] << " " << a[2] << endl;
#endif
        }

        if (scan_file_for_string(hFile, H_MARKER)) {
            cerr << " error: could not find camera parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &h[0], &h[1], &h[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "H: " << h[0] << " " << h[1] << " " << h[2] << endl;
#endif
        }

        if (scan_file_for_string(hFile, V_MARKER)) {
            cerr << " error: could not find camera parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &v[0], &v[1], &v[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "V: " << v[0] << " " << v[1] << " " << v[2] << endl;
#endif
        }
        if (scan_file_for_string(hFile, O_MARKER)) {
            gotOR = false;
        }
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &o[0], &o[1], &o[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "O: " << o[0] << " " << o[1] << " " << o[2] << endl;
#endif
        }
        if (scan_file_for_string(hFile, R_MARKER)) {
            gotOR = false;
        }
        else {
            if (fscanf(hFile, " = (%f,%f,%f", &r[0], &r[1], &r[2]) != 3) {
                cerr << " error: could not find camera parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "R: " << r[0] << " " << r[1] << " " << r[2] << endl;
#endif
        }
        fclose(hFile);
        return 0;
    }

    /////////////////////////////////////////////////////////////////////////////

    int scan_file_for_sun_angles(const char* strInFilename, 
                                 float *azimuth, float *elevation)
    {
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // first scan for the site geometry group
        if (scan_file_for_string(hFile, SITE_GEOM_MARKER)) {
            cerr << " error: could not find site derived geometry params "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
    
        // scan and grab sun azimuth and elevation
        if (scan_file_for_string(hFile, SUN_AZ_MARKER)) {
            cerr << " error: could not find sun azimuth parameters "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        } 
        else {
            if (fscanf(hFile, " = %f", azimuth) != 1) {
                cerr << " error: could not find sun azimuth parameters "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

#ifdef _DEBUG
            cerr << "Sun Azimuth: " << *azimuth << endl;
#endif
        }

        if (scan_file_for_string(hFile, SUN_EL_MARKER)) {
            cerr << " error: could not find sun elevation info "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = %f", elevation) != 1) {
                cerr << " error: could not find sun elevation info "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
#ifdef _DEBUG
            cerr << "Sun Elevation: " << *elevation << endl;
#endif
        }
        fclose(hFile);
        return 0;
    }
    /*
      int scan_file_for_solar_time(const char* strInFilename, float *solarTimeInSec)
      {

      FILE* hFile = fopen(strInFilename, "r");
      if (!hFile) {
      cerr << " error: could not open input file " 
      << strInFilename << ", exiting.\n";
      return -1;
      }
  
      // first scan for the site geometry group
      if (scan_file_for_string(hFile, LOCAL_TRUE_SOLAR_TIME )) {
      cerr << " error: could not find local true solar time "
      << "in file " << strInFilename << ", exiting.\n";
      fclose(hFile);
      return -1;
      }
 
      }
    */
    /////////////////////////////////////////////////////////////////////////////

    // Converts the 'a' vector of cahv format to pan/tilt values, assumes
    // camera out vector is the y-axis, that the coordinate system is
    // right handed, and that the z-axis points down Also, the reference
    // coordinate frame (of 'a') is assumed to be a cartesian coordinate
    // frame whose z-axis is pointing downward; this means that if tilt is
    // returned > 0, then the camera is pointing down.
    void a2pan_tilt(const double a[3], float &pan, float& tilt)
    {
        const float x = a[0];
        const float y = a[1];
        const float z = a[2];

        pan = (float)atan2((double)x, (double)y);

        // note no need to divide by norm(a) here bec. a is already normalized
        tilt = (float)atan2((double)-z, sqrt((double)(x*x+y*y)));

#ifdef _DEBUG
        cerr << "pan=" << pan << "; tilt=" << tilt << endl;
#endif
    }



    int scan_file_for_clock_start(const char* strInFilename, double &start_time)
    {
        char temp;

        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // scan for the start time marker and grab start time
        if (scan_file_for_string(hFile, START_TIME_MARKER)) {
            cerr << " error: could not find start time "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = %c%lf", &temp, &start_time) != 2) {
                cerr << " error: could not find start time "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

#ifdef _DEBUG
            cerr << "Start time: " << start_time << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }

    int scan_file_for_left_or_right(const char* strInFilename, int &left_right_id)
    {
        char temp;
        char letter;

        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // scan for the start time marker and grab start time
        if (scan_file_for_string(hFile, LR_MARKER)) {
            cerr << " error: could not find if the fame was left or right "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = %c%c", &temp,&letter) != 2) {
                cerr << " error: could not find if the fame was left or right "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

            if (letter=='L')
                left_right_id=LEFT_ID;
            else if (letter=='R')
                left_right_id=RIGHT_ID;
            else
                left_right_id=NEITHER_ID;

#ifdef _DEBUG
            cerr << "Frame id letter read: " << letter << endl;
            cerr << "left_right_id: " << left_right_id << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }

    //what is PMA?
    int scan_file_for_rover_PMA(const char* strInFilename, int &PMA)
    {
        int temp1, temp2, temp3;
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // scan for the rover motion marker and grab PMA
        if (scan_file_for_string(hFile, RMC_MARKER)) {
            cerr << " error: could not find the robot motion counter "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%d,%d,%d,%d", &temp1, &temp2, &temp3, &PMA) != 4) {
                cerr << " error: could not find the PMA "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
            //cout<<"temp1="<<temp1<<", temp2="<<temp2<<", temp3="<<temp3<<endl;
#ifdef _DEBUG
            cerr << "PMA Id in rover motion counter: " << PMA << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }

    int scan_file_for_rover_motion_counter(const char* strInFilename, int &site, int &drive)
    {
        int temp, PMA;
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // scan for the rover motion marker and grab PMA
        if (scan_file_for_string(hFile, RMC_MARKER)) {
            cerr << " error: could not find the robot motion counter "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%d,%d,%d,%d", &site, &drive, &temp, &PMA) != 4) {
                cerr << " error: could not find the PMA "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }
            //cout<<"temp1="<<temp1<<", temp2="<<temp2<<", temp3="<<temp3<<endl;
#ifdef _DEBUG
            cerr << "PMA Id in rover motion counter: " << PMA << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }

    int scan_file_for_origin_offset(const char* strInFilename, double offset[3])
    {
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // first scan for the rover coordinate system group
        if (scan_file_for_string(hFile, ROVER_COORDINATE_MARKER)) {
            cerr << " error: could not find the rover coordinates "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }


        // scan for the offset marker and grab offset vector
        if (scan_file_for_string(hFile, OFFSET_MARKER)) {
            cerr << " error: could not find the origin offset "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%lf,%lf,%lf", offset, offset+1, offset+2) != 3) {
                cerr << " error: could not find the offset "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

#ifdef _DEBUG
            cerr << "Origin offset vector: (" << offset[0] << ", " << offset[1] << ", " << offset[2] << ")" << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }


    int scan_file_for_origin_rotation(const char* strInFilename, double quaternion[4])
    {
        FILE* hFile = fopen(strInFilename, "r");
        if (!hFile) {
            cerr << " error: could not open input file " 
                 << strInFilename << ", exiting.\n";
            return -1;
        }

        // first scan for the rover coordinate system group
        if (scan_file_for_string(hFile, ROVER_COORDINATE_MARKER)) {
            cerr << " error: could not find the rover coordinates "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }

        // scan for the offset marker and grab rotation
        if (scan_file_for_string(hFile, ROTATION_MARKER)) {
            cerr << " error: could not find the origin rotation "
                 << "in file " << strInFilename << ", exiting.\n";
            fclose(hFile);
            return -1;
        }
        else {
            if (fscanf(hFile, " = (%lf,%lf,%lf,%lf", quaternion, quaternion+1, quaternion+2, quaternion+3) != 4) {
                cerr << " error: could not find the origin rotation "
                     << "in file " << strInFilename << ", exiting.\n";
                fclose(hFile);
                return -1;
            }

#ifdef _DEBUG
            cerr << "Origin offset vector: (" << quaternion[0] << ", " << quaternion[1] << ", " << quaternion[2] 
                 << ", " << quaternian[3] << ")" << endl;
#endif
        }

        fclose(hFile);
        return 0;
    }


    int scan_file_for_image_size( const char* strInFilename, long &nCols, long &nRows )
    {

        // parse the pds file to get the top level node odltree
        ODLTREE odltree = OaParseLabelFile((char *)strInFilename, 
                                           (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);

        if (odltree == NULL) { 
            OaReportError( (char*) "Error in OaParseLabelFile()" );
            return -1; 
        }

        // grab the image node from the top level node
        ODLTREE image_node = OdlFindObjDesc(odltree, (char*) "IMAGE", NULL, NULL,
                                            ODL_RECURSIVE_DOWN, ODL_TO_END);

        if (image_node == NULL) {
            OaReportError( (char*) "Error in OaParseLabelFile()" );
            return -1;
        }

        nRows = get_pds_label_value(image_node, "LINES");
        nCols = get_pds_label_value(image_node, "LINE_SAMPLES");


        OdlFreeTree( odltree );

        return 0;

    }
    /////////////////////////////////////////////////////////////////////////////
    int getExposureDuration(const char* strInFilename, double *exposureDuration)
    {

        // parse the pds file to get the top level node odltree
        ODLTREE odltree = OaParseLabelFile((char *)strInFilename, 
                                           (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);

        if (odltree == NULL) { 
            OaReportError( (char*) "Error in OaParseLabelFile()" );
            return -1; 
        }
  
        KEYWORD *groupKeyword, *parameterKeyword;
        groupKeyword = OdlFindKwd(odltree, "GROUP", "INSTRUMENT_STATE_PARMS", 1, 
                                  ODL_THIS_OBJECT);
        if (groupKeyword == NULL) {
            cerr << "ERROR in getExposureDuration(): Could not find INSTRUMENT_STATE_PARMS!"
                 << endl;
            return -1;
        }
  
        parameterKeyword = OdlNextKwd(groupKeyword, "EXPOSURE_DURATION", "*",
                                      1, ODL_THIS_OBJECT);
        if (parameterKeyword == NULL) {
            cerr << "ERROR in getExposureDuration(): Could not find EXPOSURE_DURATION"
                 << endl;
            return -1;
        }

        char *conversionPos = 0;
        char *keywordValue = OdlGetKwdValue(parameterKeyword);

        if (strcmp(keywordValue, "\"NULL\"") == 0)
            return -1;

        *exposureDuration = strtod(keywordValue, &conversionPos);
        if (conversionPos == keywordValue)
        {
            cerr << "ERROR in getExposureDuration(): invalid EXPOSURE_DURATION value"
                 << endl;
            return -1;
        }

        return 0;
    }
    //reads the entire data buffer in unsigned char and returns number of bands, bits and image size.
    int PDSReadImage(string strInFilename, unsigned char** r_pDataBuffer, 
                     long *r_nCols, long *r_nRows, long *r_nBands, long *r_nBits)
    {
	
        // parse the pds file to get the top level node odltree
        ODLTREE odltree = OaParseLabelFile((char *)strInFilename.c_str(), 
                                           (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);

        if (odltree == NULL) { 
            OaReportError((char*) "Error in OaParseLabelFile()");
            return -1; 
        }
 
        // grab the image node from the top level node
        ODLTREE image_node = OdlFindObjDesc(odltree, (char*) "IMAGE", NULL, NULL,
                                            ODL_RECURSIVE_DOWN, ODL_TO_END);

        if (image_node == NULL) {
            OaReportError( (char*) "Error in OaParseLabelFile()");
            return -1;
        }

        ///////////////////////////////////////////////////////////////////////////
        // process label part of PDS file
        ///////////////////////////////////////////////////////////////////////////

        const long nBands = get_pds_label_value(image_node, "BANDS");
        *r_nBands = nBands;
        /*
          if (nBands != 1) {
          cerr << "read_pds_image: only single-band images are supported, exiting.\n";
          return -1;
          }
        */
        const long nBits =  get_pds_label_value(image_node, "SAMPLE_BITS");
        *r_nBits = nBits;
        /* 
           if (nBits != INPUT_BITS_PER_PIXEL) {
           cerr << "read_pds_image: unsupported bits/pixel, exiting." << endl;
           return -1;
           }
        */
        const long nRows = get_pds_label_value(image_node, "LINES");
        const long nCols = get_pds_label_value(image_node, "LINE_SAMPLES");
        *r_nCols = nCols;
        *r_nRows = nRows;
        printf("%ld %ld\n", nRows, nCols);
   
        ///////////////////////////////////////////////////////////////////////////
        // process data part of PDS file (image data)
        ///////////////////////////////////////////////////////////////////////////

        unsigned char *pDataBuffer = new unsigned char[nRows*nCols*nBands*nBits/8];
 
        for(int band=0; band<nBands; ++band) {

            cout<<"band="<<band<<endl;
            // open the image (the '1' argument is the band #, assuming band #'s
            // are 1-indexed, not zero-indexed.. -DT)
            OA_OBJECT image_handle = OaOpenImage(image_node, 1);

            if (image_handle == NULL) {
                OaReportError((char*) "Error in OaOpenImage()");
                return -1;
            }

            // NOTE: batch version -- no copying done, you may want to read the
            // image by chunks here, if it's huge, using OaReadPartialImage()
            OA_OBJECT image_obj = OaReadImage(image_node, band+1);
            if (image_obj == NULL) {
                OaReportError((char*) "Error in OaReadImage()");
                return -1;
            }
      
            switch (nBits){
            case 8: 
                {
                    for (long i = 0; i < nRows*nCols; ++i){
                        pDataBuffer[band*nRows*nCols+i]=image_obj->data_ptr[i];
                    }
                    break;
                }
            case 16:
                {

                    for (long i = 0; i < nRows*nCols*2; ++i){
                        pDataBuffer[band*nRows*nCols*2+i]=image_obj->data_ptr[i];
                    }
                    break;
                }
            case 32:
                {
                    for(long i = 0; i < nRows*nCols*4; ++i){
                        pDataBuffer[band*nRows*nCols*4+i]=image_obj->data_ptr[i];
                    }
                    break;
                }
            default:
                {
                    OaReportError((char*) "Error in PDSReadImage(): nBits is not 8/16/32"); 
                }    
            }

            OaCloseImage(image_handle);
            OaDeleteObject(image_obj);

        }
        *r_pDataBuffer = pDataBuffer;
        return 0;
    }

    void Convert16BitTo8BitImage(unsigned short int* inputBuffer, unsigned char *outputBuffer, int nCols, int nRows)
    {
 
        //determine the min max value of the inputBuffer
        long minVal = 10000000;
        long maxVal = -10000000; 
        long currVal;
        for (int i = 0; i <nRows; i++){
            for (int j = 0; j <nCols; j++){
                currVal = inputBuffer[i*nCols+j];
     
                if (currVal<minVal){
                    minVal = currVal;
                }
                if (currVal>=maxVal){
                    maxVal = currVal;
                }
            }
        }
        //printf("minVal = %ld, maxVal = %ld\n", minVal, maxVal);

        //determine the scalling and offset 
        float scalingFactor = 255.0/maxVal;
   
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                float val = scalingFactor*inputBuffer[i*nCols+j];
                if (val > 255){
                    val = 255;
                }
                if (val < 0){
                    val = 0;
                }
                outputBuffer[i*nCols+j] = (unsigned char)val; 
            }
        } 
   
    }
    uint8_t *
    create_8bit_image_from_16bit_image(const uint16_t *valueBuffer,
                                       const uint32_t numValues)
    {
        const static uint16_t kMax16BitValue = numeric_limits<uint16_t>::max();
        const static uint8_t kMax8BitValue = numeric_limits<uint8_t>::max();
  
        uint8_t *destBuffer = new uint8_t[numValues];
        if (destBuffer == 0)
        {
            cout << "FATAL ERROR in create_8bit_image_from_16bit_image():"
                 << " couldn't allocate destBuffer. Aborting." << endl;
            return 0;
        }

        double scale = double(kMax8BitValue) / double(kMax16BitValue);
        cout<<"scale="<<scale<<endl;
        for (uint32_t i = 0; i < numValues; i++)
        {
            //cout<<"vBuffer="<<((uint16_t*)valueBuffer)[i]<<endl;
            float value = float(valueBuffer[i]) * scale;
            /*
              cout<<"value="<<value<<", "<<round(value)<<endl;
              if (value>kMax8BitValue){
              destBuffer[i]=kMax8BitValue;
              }
              else{
              destBuffer[i]=uint8_t(round(value));
              cout<<"dBuffer="<<destBuffer[i]<<endl;
              }
            */
            destBuffer[i] = ((value > kMax8BitValue) ?
                             kMax8BitValue : /*uint8_t*/static_cast<uint8_t>(round(value)));
            //cout<<"dBuffer="<<uint8_t(destBuffer[i])<<endl;
        }

        return destBuffer;
    }

    void WriteCAHVORModel(string strCalFilename, float c[3], float a[3], float h[3], float v[3], float o[3], float r[3])
    {
        // write the CAHV params to a camera calibration file for stereo pipeline
        printf("CalFilename = %s\n", strCalFilename.c_str());
        FILE *cameraCalFile = fopen(strCalFilename.c_str(), "w");
  
        fprintf(cameraCalFile, "C = %14.6f %14.6f %14.6f\n", c[0], c[1], c[2]);
        fprintf(cameraCalFile, "A = %14.6f %14.6f %14.6f\n", a[0], a[1], a[2]);
        fprintf(cameraCalFile, "H = %14.6f %14.6f %14.6f\n", h[0], h[1], h[2]);
        fprintf(cameraCalFile, "V = %14.6f %14.6f %14.6f\n", v[0], v[1], v[2]);
        fprintf(cameraCalFile, "O = %14.6f %14.6f %14.6f\n", o[0], o[1], o[2]);
        fprintf(cameraCalFile, "R = %14.6f %14.6f %14.6f\n", r[0], r[1], r[2]);
        fprintf(cameraCalFile, "W = 1024\n");
        fprintf(cameraCalFile, "H = 1024\n");

        fclose(cameraCalFile);
    }

    void WriteCAHVModel(string strCalFilename, float c[3], float a[3], float h[3], float v[3])
    {
        // write the CAHV params to a camera calibration file for stereo pipeline
        printf("CalFilename = %s\n", strCalFilename.c_str());
        FILE *cameraCalFile = fopen(strCalFilename.c_str(), "w");
  
        fprintf(cameraCalFile, "C = %14.6f %14.6f %14.6f\n", c[0], c[1], c[2]);
        fprintf(cameraCalFile, "A = %14.6f %14.6f %14.6f\n", a[0], a[1], a[2]);
        fprintf(cameraCalFile, "H = %14.6f %14.6f %14.6f\n", h[0], h[1], h[2]);
        fprintf(cameraCalFile, "V = %14.6f %14.6f %14.6f\n", v[0], v[1], v[2]);
        fprintf(cameraCalFile, "W = 1024\n");
        fprintf(cameraCalFile, "H = 1024\n");
        /*
          W =    1024
          H =    1024
        */
        fclose(cameraCalFile);
    }

  #if 0
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
    
    unsigned char* read_pgm(string filename, long *numCols, long *numRows)
    {

      FILE* file = fopen(filename.c_str(), "r");
      if (!file) return(false);
      
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
  #endif
}
