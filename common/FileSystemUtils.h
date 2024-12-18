#ifndef FILESYSTEM_UTIL_H_
#define FILESYSTEM_UTIL_H_

#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>

//
#include <vector>
#include <string>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <limits>
#include <unistd.h>
#include <stdlib.h>

namespace at
{
  inline bool MakeDirectory(std::string dirPath)
  {
    if ((mkdir(dirPath.c_str(), S_IRWXU|S_IRWXG|S_IRWXO) == -1) &&
	(errno != EEXIST))
    {
      std::string command = "mkdir(\"" + dirPath + "\", S_IRWXU|S_IRWXG|S_IRWXO)";
      perror(command.c_str());
      return false;
    }
    else if (errno == EEXIST)
    {
      struct stat stat_buffer;
      // We don't error check here since we know file exists
      stat(dirPath.c_str(), &stat_buffer);
      if (!S_ISDIR(stat_buffer.st_mode))
      {
	std::cerr << "ERROR in at::CreateDirectory(): " << dirPath
		  << " exists and is not a directory." << std::endl;
	return false;
      }
    }
    return true;
  }
   

  bool is_file_exist(const char *fileName);

  //this is a cool function that deserves to be moved to common
  //returns in variable files the filenames in directory dir
  // This function will only return files with IMG extensions
  // Names are returned in an arbitrary order, sorting is recommended
  int getFileNames (std::string dir, std::vector<std::string> &files);

  // This function will return files with all extensions
  // includeDirs will include directory names from results if set to true
  // Names are returned in an arbitrary order, sorting is recommended
  int getFileNames (std::string dir, std::vector<std::string> &files, bool includeDirs);
  int GetImageFileNames (std::string dir, std::vector<std::string> &files);
}

#endif	// FILESYSTEM_UTIL_H_
