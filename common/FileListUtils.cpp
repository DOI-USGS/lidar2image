// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>		// for atof()

#include "FileListUtils.h"
#include "StringUtils.h"

using namespace std;

namespace at
{

     void PrintOverlapList( const std::vector<int>& overlapIndices ) {
        cout << "numOverlapping = " << overlapIndices.size() << endl;
        std::vector<int>::const_iterator i;
        for( i = overlapIndices.begin(); i != overlapIndices.end(); ++i){
            cout << "index: " << *i << endl;
        }
    }

    void SaveOverlapList( const string& filename, const std::vector<int>& overlapIndices ) {
        WriteVectorTo( filename, overlapIndices );
    }

    void SaveOverlapList( const string& filename, const std::vector<std::string>& filenames ) {
        WriteVectorTo( filename, filenames );
    }

    int ReadOverlapList( const string& filename, std::vector<int>& overlapIndices ) {
        overlapIndices = ReadVectorFrom<int>( filename );
        return 1;
    }

    void ReadFileList(const string& filePath, vector<string>& fileArray )
    {
      string line;
      ifstream file (filePath.c_str());
      fileArray.resize(0);
      if (file.is_open()){ 
	while (!file.eof()){
	  std::getline(file, line);
	  if (!line.empty())
	    fileArray.push_back(line);
	}
	file.close();
      }      
    }
  
    // Read list of double from textfile
    vector<double> ReadFileListDouble(string filename)
    {
      vector<double> timeStampVector;
      timeStampVector.resize(0);
      string line;
      ifstream file (filename.c_str());
      
      if (file.is_open()){ 
	while (!file.eof()){
	  std::getline(file, line);
	  if (!line.empty())
	    timeStampVector.push_back(atof(line.c_str()));
	}
	file.close();
      }
      return timeStampVector;
   }
  
   vector<string> AccessDataFilesFromInput(const string& inputFile) 
   {
      const vector<string> v(1,inputFile);
      return AccessDataFilesFromInput( v );
   }
  /*
    vector<string> AccessDataFilesFromInput(const vector<string>& inputFiles) 
    {
      vector<string> v;
      if( inputFiles.empty() ){ 
	cout<< "There are no input files to read."<<endl;
      }
      for( vector<string>::const_iterator i = inputFiles.begin();
	   i != inputFiles.end();
	   ++i ){
	fs::path p( *i );
	if( p.extension() == string(".txt") ){
	  vector<string> list = ReadVectorFrom<string>( p );
	  v.insert( v.end(), list.begin(), list.end() );
	}
	else{ v.push_back( *i ); }
      }
      return v;
    }
  */
  
  vector<string> AccessDataFilesFromInput(const vector<string>& inputFiles) 
  {
    vector<string> v;
    if( inputFiles.empty() ){ 
      cout<< "There are no input files to read."<<endl;
    }
    for (unsigned int i=0; i < inputFiles.size(); i++){
      string extension = GetFilenameExt(inputFiles[i]);
      if (extension.compare("txt")==0){
	vector<string> thisFilenameArray;
	ReadFileList(inputFiles[i], thisFilenameArray);
	for (unsigned int j=0; j <thisFilenameArray.size(); j++){
	  v.push_back(thisFilenameArray[j]);
	}
      }
      else{
	v.push_back(inputFiles[i]);
      }
    }
    return v;
  }
  
  vector<string> AccessDataFilesFromInput(const string &inputFilename, string filenameExtension) 
  { 
    //read the mastcam surface filename
    vector<string >filenameArray;
    string extension = GetFilenameExt(inputFilename);
    cout << "input filename = " << inputFilename << ", "
	 << "extension = " << extension << endl;
    
    if (extension.compare(filenameExtension) == 0) {
      filenameArray.push_back(inputFilename);
    }  
    else if (extension.compare("txt") == 0) {
      ReadFileList(inputFilename, filenameArray);
    }
  
    return filenameArray;
  }

  // Send complete file extensions. Include '.'. Ex: ".txt"
  vector<string> AccessDataFilesFromInput(const string &inputFilename, const string &filenameExtension,
					  const string &listnameExtension) 
  {

    vector<string >filenameArray;

    // Check for list file or single file
    std::size_t found = inputFilename.find(listnameExtension);

    if (found!=std::string::npos)
    // This is a list of files
    {
      ReadFileList(inputFilename, filenameArray);
    }  
    else 
    {
      // Check if it is a single file
      found = inputFilename.find(filenameExtension);
      if (found!=std::string::npos) {
      filenameArray.push_back(inputFilename);
      }
    }
    return filenameArray;
  }
 
}
