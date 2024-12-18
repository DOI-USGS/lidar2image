#ifndef FILE_LIST_UTILS_H
#define FILE_LIST_UTILS_H

// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
//

#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

namespace at
{

    void ReadFileList( const std::string&, std::vector<std::string>& );
    std::vector<double> ReadFileListDouble(std::string filename);
    void PrintOverlapList( const std::vector<int>& );
    void SaveOverlapList( const std::string&, const std::vector<int>& );
    void SaveOverlapList( const std::string&, const std::vector<std::string>& );
    int  ReadOverlapList( const std::string&, std::vector<int>&);
    
    
    //determine if inputFiles is a text file containing a list of input files, one input file
    //or a set of input files
    //by reading the inputFiles extension. A text file containing a filename list *must*
    //have extension .txt 
    std::vector<std::string> AccessDataFilesFromInput(const std::string& );
    std::vector<std::string> AccessDataFilesFromInput(const std::vector<std::string>& );
    std::vector<std::string> AccessDataFilesFromInput(const std::string& inputFilename,
						      std::string filenameExtension); 

    // Used to specify list and file extensions. Useful when the normal file ends in .txt
    // Example usage: 
    //normal file name = test.txt
    //filenameExtension = .txt
    //list file name = testAbunch_list.txt
    //listnameExtension = _list.txt
    // Send complete file extensions. Include '.'. Ex: ".txt"
    std::vector<std::string> AccessDataFilesFromInput(const std::string &inputFilename,
						      const std::string &filenameExtension,
						      const std::string &listnameExtension);
    
    /* Stream manipulator to ignore until the end of line
     * Use like: std::cin >> ignoreLine;
     */
    template <class charT, class traits>
    inline
    std::basic_istream<charT, traits>&
    ignoreLine (std::basic_istream<charT,traits>& stream)
    {
        // skip until end of line
        stream.ignore( std::numeric_limits<int>::max(), stream.widen('\n') );
        return stream;
    }


    /* Stream manipulator to ignore one character
     * Use like: std::cin >> ignoreOne;
     */
    template <class charT, class traits>
    inline
    std::basic_istream<charT, traits>&
    ignoreOne(std::basic_istream<charT,traits>& stream)
	{
        stream.ignore();

        return stream;
	}

    template <class T>
    inline std::vector<T> ReadVectorFrom( std::istream& stream ){
        std::vector<T> v;
        T item;
        while( stream >> item ){ v.push_back( item ); }
        return v;
    }

    template <class T>
    inline std::vector<T> ReadVectorFrom( const std::string& s ){
      std::ifstream file( s.c_str() );
      if (!file){
          std::cout<<"Can't open file " << s <<std::endl;
      }
      return ReadVectorFrom<T>( file );
    }
      
    template <class T>
    inline std::ostream& WriteVectorTo( std::ostream& stream, const std::vector<T>& v){
        if( v.empty() ){ return stream; }
        for( typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i ){
            stream << *i << std::endl;
        }
        return stream;
    }

    template <class T>
    inline void WriteVectorTo( const std::string& s, const std::vector<T>& v){
        if( v.empty() ){
            std::cout<< "WriteVectorTo: the input vector is empty."<<std::endl;
        }

        std::ofstream file( s.c_str() );
        if( !file ) {
            std::cout<<"Can't open file \"" << s << "\""<<std::endl;
        }

        WriteVectorTo( file, v );
        file.close();
        return;
    }

  
}

#endif	// FILE_LIST_UTILS_H

