// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "StringUtils.h"

using namespace std;

namespace at
{

    // From: http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
    // These functions split a string given a delimiter
    // Should these go to common/string_util? The current FindAndSplit doesn't get the last item on the line
    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
    {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
	elems.push_back(item);
      }
      return elems;
    };
     
    std::vector<std::string> split(const std::string &s, char delim/* = ' '*/)
    {
      std::vector<std::string> elems;
      split(s, delim, elems);
      return elems;
    };
  
    vector<std::string> FindAndSplit(std::string& tInput, std::string tFind)
    {
        vector<std::string> elems;
        size_t uFindLen = tFind.length(); 
  
        if( uFindLen != 0 ){

            int index = tInput.find(tFind);

            while (index != -1){
                string elem(tInput, 0, index);
		cout<<"elem="<<elem<<endl;
                tInput.erase(0, index+uFindLen);
                elems.push_back(elem);
                index = tInput.find(tFind); 
            }
        }
        return elems;
    }

    void FindAndReplace( std::string& tInput, std::string tFind, std::string tReplace ) 
    { 

        size_t uPos = 0; 
        size_t uFindLen = tFind.length(); 
        size_t uReplaceLen = tReplace.length();
  
        if( uFindLen != 0 ){
            for( ;(uPos = tInput.find( tFind, uPos )) != std::string::npos; ){
                tInput.replace( uPos, uFindLen, tReplace );
                uPos += uReplaceLen;
            }
        }

    }

    ///Returns the file extension
    std::string GetFilenameExt(std::string const& filename) {
        std::string result = filename;
        int index = result.rfind(".");
        if (index != -1)
            result.erase(0, index+1);
        return result;
    }

    /// Erases a file suffix if one exists and returns the base string
    std::string GetFilenameNoExt(std::string const& filename) {
        std::string result = filename;
        int index = result.rfind(".");
        if (index != -1)
            result.erase(index, result.size());
        return result;
    }

    // Erases a file path if one exists and returns the base string 
    std::string GetFilenameNoPath(std::string const& filename) {
        std::string result = filename;
        int index = result.rfind("/");
        if (index != -1)
            result.erase(0, index+1);
        return result;
    }

    std::string GetDirname(std::string const& filename)
    {
        std::string dirname=filename;
        int index = dirname.rfind("/");
        if (index != -1){
	  dirname.erase(index/*+1*/, filename.size());
	}
        else{
	  dirname="";
        }
        return dirname;
    }

    //returns the string within original string and in between string before and string after
    string findString(string origString, string before, string after)
    {
        string resultString;
        size_t beforeIndex = origString.rfind(before)+before.size();
        size_t afterIndex  = origString.rfind(after);
        int resultLength = afterIndex-beforeIndex;
        resultString.resize(resultLength);
        //cout<<"beforeIndex="<<beforeIndex<<", afterIndex="<<afterIndex<<", length="<<resultLength<<", origLength="<<origString.size()<<endl;
        origString.copy((char*)resultString.c_str(), resultLength, beforeIndex);
        return resultString;
    }

    std::string ZeroPadNumber(int num, int stringWidth)
    {
      std::ostringstream ss;
      ss << std::setw( stringWidth) << std::setfill( '0' ) << num;
      return ss.str();
    }


    std::string NumberToString(int number)
    {
      string str;
      std::ostringstream ss;
      ss<<number;
      return ss.str();
    }
  
}
