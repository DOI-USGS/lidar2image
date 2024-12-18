#ifndef TILING_H
#define TILING_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

namespace at
{
  
    struct BoundingBox
    {
      std::vector<float> first;
      std::vector<float> last;
    };  

    //these is the return structure
    struct TilingParams
    {
      //pixel coordinate in the background image
      //corresponding to the tiles
      int back_xl; 
      int back_yt;
      int back_xr;
      int back_yb;
      
      //pixel coordinates in the background image
      //corresponding to area of overlap with 
      //foreground images 
      double fore_xl;
      double fore_yt;
      double fore_xr;
      double fore_yb;
      
      int horTileIndex;
      int verTileIndex;
    };
    
    struct bbox
    {
      double xl;
      double xr;
      double yt;
      double yb;
    };
    
    struct whBox
    {
      int xl;
      int yt;
      int width;
      int height;
    };

    struct TileConfig
    {
      int   imageROILeft;
      int   imageROITop;
      int   imageROIRight;
      int   imageROIBottom;
      int   tileWidth;
      int   tileHeight;
      int   overlapHor;
      int   overlapVer;
      int   paddingTop;
      int   paddingLeft;
      int   paddingRight;
      int   paddingBottom;
      float pyrResampleFactor;
    };

    class Tiling
    {
    public:
      //Constructor
      Tiling();
      Tiling(std::string configFilename);
      Tiling(struct TileConfig tileConfigParams);       
      
      //Destructor
      ~Tiling();
      
      //this member variable will be replaced with the variables below
      TileConfig m_tileParams;
      
      //these are the class variable that will be used to replace the structure above 
      int   m_imageROILeft;
      int   m_imageROITop;
      int   m_imageROIRight;
      int   m_imageROIBottom;
      int   m_tileWidth;
      int   m_tileHeight;
      int   m_overlapHor;
      int   m_overlapVer;
      int   m_paddingTop;
      int   m_paddingLeft;
      int   m_paddingRight;
      int   m_paddingBottom;
      float m_pyrResampleFactor;
      
      //this is old style, will be replaced
      std::vector<whBox> m_tiles;
      
      //will be replaced with this
      std::vector<std::vector<struct TilingParams> > m_pyrTileParamsArray;
      std::vector<float> m_offsetPix;
      
      //Methods
      void setDefaultParameters();
      void setParams(struct TileConfig tileConfigParams);
      void printTilingParams();
      
      //these three functions will be combined into one - START
      void process(int nRows, int nCols);
      void processQuadTree(struct bbox bboxPix,int pyramidMode);

      //void processPyramid(struct bbox bboxPix);
      //void processSubPyramid(struct bbox bboxPix);
      //these three functions will be combined into one - END
      
      void saveTileBB(std::string filename);
      struct TilingParams getBBox(int l, int hor, int ver);
      
    };
    
    bool isRectangle(struct bbox bBox);
    struct bbox computeBoxIntersection(struct bbox box1, struct bbox box2);
    void printBBox3(struct BoundingBox bbox3, std::string bBoxID);

    
  
  
}

#endif	// TILING_H
