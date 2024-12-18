#ifndef STRING_UTIL_H
#define STRING_UTIL_H

// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
//

#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>      // std::setw

using namespace std;

namespace at
{
  std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);    
  std::vector<std::string> split(const std::string &s, char delim = ' ');
  void FindAndReplace( std::string& tInput, std::string tFind, std::string tReplace );
  vector<std::string> FindAndSplit(std::string& tInput, std::string tFind);
  std::string GetFilenameNoExt(std::string const& filename);
  std::string GetFilenameNoPath(std::string const& filename);
  std::string GetFilenameExt(std::string const& filename);
  std::string GetDirname(std::string const& filename);
  std::string findString(std::string origString, std::string before, std::string after);
  std::string ZeroPadNumber(int num, int stringWidth);
  std::string NumberToString(int number);
  
  template<typename T> T StringToNumber(const std::string& numberAsString)
  {
    T valor;    
    std::stringstream stream(numberAsString);
    stream >> valor;
    return valor;
  }
}

#endif	// STRING_UTIL_H
