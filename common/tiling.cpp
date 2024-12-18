#include "Tiling.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>

using namespace std;

namespace at 
{

  Tiling::Tiling()
  {
    setDefaultParameters();
  }

  Tiling::Tiling(std::string configFilename)
  {

    string line;
    ifstream fin;
    string identifier;

    setDefaultParameters();
    
    fin.open(configFilename.c_str());

    if(fin.is_open()){

      cout<< "Tiling configuration file "<<configFilename<<" found."<<endl;
      while (!fin.eof()){
	std::getline(fin, line);
	if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines   
 
	  istringstream lineStream(line);
	  string keyword("");
	  lineStream >> keyword;
	   
	  if (keyword.compare(string("IMAGE_ROI_LEFT"))==0){	 
	    lineStream >>m_tileParams.imageROILeft;
	  }
	  if (keyword.compare(string("IMAGE_ROI_TOP"))==0){	 
	    lineStream >>m_tileParams.imageROITop;
	  }	   
	  if (keyword.compare(string("IMAGE_ROI_RIGHT"))==0){	 
	    lineStream >>m_tileParams.imageROIRight;
	  }
	  if (keyword.compare(string("IMAGE_ROI_BOTTOM"))==0){	 
	    lineStream >>m_tileParams.imageROIBottom;
	  }
	  if (keyword.compare(string("TILE_WIDTH"))==0){	 
	    lineStream >>m_tileParams.tileWidth;
	  }
	  if (keyword.compare(string("TILE_HEIGHT"))==0){	 
	    lineStream >>m_tileParams.tileHeight;
	  }
	  if (keyword.compare(string("TILE_OVERLAP_X"))==0){	 
	    lineStream >>m_tileParams.overlapHor;
	  }
	  if (keyword.compare(string("TILE_OVERLAP_Y"))==0){	 
	    lineStream >>m_tileParams.overlapVer;
	  }
	  if (keyword.compare(string("TILE_PADDING_TOP"))==0){	 
	    lineStream >>m_tileParams.paddingTop;
	  }
	  if (keyword.compare(string("TILE_PADDING_LEFT"))==0){	 
	    lineStream >>m_tileParams.paddingLeft;
	  }
	  if (keyword.compare(string("TILE_PADDING_RIGHT"))==0){	 
	    lineStream >>m_tileParams.paddingRight;
	  }
	  if (keyword.compare(string("TILE_PADDING_BOTTOM"))==0){	 
	    lineStream >>m_tileParams.paddingBottom;
	  }
	  if (keyword.compare(string("PYRAMID_RESAMPLE_FACTOR"))==0){	 
	    lineStream >>m_tileParams.pyrResampleFactor;
	  }

	}
      }

    }
    else{
      cout<<"Warning! Configuration file not found. Using default values."<<endl;
    }

    printTilingParams();
  }

  Tiling::Tiling(struct TileConfig tileConfigParams)
  {
    m_tileParams.imageROILeft = tileConfigParams.imageROILeft;
    m_tileParams.imageROITop = tileConfigParams.imageROITop;
    m_tileParams.imageROIRight = tileConfigParams.imageROIRight;
    m_tileParams.imageROIBottom = tileConfigParams.imageROIBottom;
    m_tileParams.tileWidth = tileConfigParams.tileWidth;
    m_tileParams.tileHeight = tileConfigParams.tileHeight;
    m_tileParams.overlapHor = tileConfigParams.overlapHor;
    m_tileParams.overlapVer = tileConfigParams.overlapVer;
    m_tileParams.paddingTop = tileConfigParams.paddingTop;
    m_tileParams.paddingLeft = tileConfigParams.paddingLeft;
    m_tileParams.paddingRight = tileConfigParams.paddingRight;
    m_tileParams.paddingBottom = tileConfigParams.paddingBottom;
    m_tileParams.pyrResampleFactor = tileConfigParams.pyrResampleFactor;
  }
   
  Tiling::~Tiling()
  {
    m_tiles.resize(0);
  }
  
  void Tiling::setDefaultParameters()
  {
    m_tileParams.imageROILeft = 0;
    m_tileParams.imageROITop = 0;
    m_tileParams.imageROIRight = 0;
    m_tileParams.imageROIBottom = 0;
    m_tileParams.tileWidth = 512;
    m_tileParams.tileHeight = 512;
    m_tileParams.overlapHor = 0;
    m_tileParams.overlapVer = 0;
    m_tileParams.paddingTop = 0;
    m_tileParams.paddingLeft = 0;
    m_tileParams.paddingRight = 0;
    m_tileParams.paddingBottom = 0;
    m_tileParams.pyrResampleFactor = 1;
  }

  void Tiling::setParams(struct TileConfig tileConfigParams)
  {
    m_tileParams.imageROILeft = tileConfigParams.imageROILeft;
    m_tileParams.imageROITop = tileConfigParams.imageROITop;
    m_tileParams.imageROIRight = tileConfigParams.imageROIRight;
    m_tileParams.imageROIBottom = tileConfigParams.imageROIBottom;
    m_tileParams.tileWidth = tileConfigParams.tileWidth;
    m_tileParams.tileHeight = tileConfigParams.tileHeight;
    m_tileParams.overlapHor = tileConfigParams.overlapHor;
    m_tileParams.overlapVer = tileConfigParams.overlapVer;
    m_tileParams.paddingTop = tileConfigParams.paddingTop;
    m_tileParams.paddingLeft = tileConfigParams.paddingLeft;
    m_tileParams.paddingRight = tileConfigParams.paddingRight;
    m_tileParams.paddingBottom = tileConfigParams.paddingBottom;
    m_tileParams.pyrResampleFactor = tileConfigParams.pyrResampleFactor;
  }

  void Tiling::printTilingParams()
  {
    cout<<"IMAGE_ROI_LEFT "<<m_tileParams.imageROILeft<<endl;
    cout<<"IMAGE_ROI_TOP "<<m_tileParams.imageROITop<<endl;
    cout<<"IMAGE_ROI_RIGHT "<<m_tileParams.imageROIRight<<endl;
    cout<<"IMAGE_ROI_BOTTOM "<<m_tileParams.imageROIBottom<<endl;
    cout<<"TILE_WIDTH "<<m_tileParams.tileWidth<<endl;
    cout<<"TILE_HEIGHT "<<m_tileParams.tileHeight<<endl;
    cout<<"TILE_OVERLAP_X "<<m_tileParams.overlapHor<<endl;
    cout<<"TILE_OVERLAP_Y "<<m_tileParams.overlapVer<<endl;
    cout<<"TILE_PADDING_TOP "<<m_tileParams.paddingTop<<endl;	 
    cout<<"TILE_PADDING_LEFT "<<m_tileParams.paddingLeft<<endl;
    cout<<"TILE_PADDING_RIGHT "<<m_tileParams.paddingRight<<endl;	 
    cout<<"TILE_PADDING_BOTTOM "<<m_tileParams.paddingBottom<<endl; 
    cout<<"PYRAMID_RESAMPLE_FACTOR "<<m_tileParams.pyrResampleFactor<<endl;
  }
 
  void Tiling::saveTileBB(string filename)
  {
    ofstream fout;
    fout.open(filename.c_str());
    for(unsigned int i=0; i<m_tiles.size(); i++){
      fout << i << "\t" << m_tiles[i].xl << " " << m_tiles[i].yt << " " << m_tiles[i].width << " " << m_tiles[i].height << endl;
    }
    fout.close();
  }
  
  //creates a set of overlapping tiles
  //tiles have no padding, and no pyramid is generated
  void Tiling::process(int nRows, int nCols)
  {
    whBox roi;
    int r,c,x,y,oldC,oldR;
    int xOver, yOver;
    
    x = m_tileParams.tileWidth;
    y = m_tileParams.tileHeight;
    xOver = m_tileParams.overlapHor;
    yOver = m_tileParams.overlapVer;
        
    m_tiles.clear();
    
    //Check if Tile too large
    if( x > nRows || x < 0){
      x = nCols;
    }
    if( y > nRows || y < 0 ){
      y = nRows;
    }
    
    for (r = 0; r < nRows; r += (y-yOver)){   
      for (c = 0; c < nCols; c += (x-xOver)){
	
	roi.width = x;
	roi.height = y;

	if(c+x > nCols){
	  roi.xl = nCols-x;
	}
	else{
	  roi.xl = c;
	}
	if(r+y > nRows){
	  roi.yt = nRows-y;
	}
	else{
	  roi.yt = r;
	}
	
	m_tiles.push_back(roi);
	
	if( c+x >= nCols){
	  break;
	}
	
      }
      
      if( r+y >= nRows ){
	break;
      }
      
    }
  }

  void Tiling::processQuadTree(struct bbox bBoxPixel, int pyramidMode)
  {

    int imgWidth = round(bBoxPixel.xr-bBoxPixel.xl);
    int imgHeight = round(bBoxPixel.yb-bBoxPixel.yt);
    cout<<"Tiling.cc:: processQuadTree(): imgWidth="<<imgWidth<<", imgHeight="<<imgHeight<<endl;
    //the effective tile size without padding
    int ROIImgWidth =  imgWidth  - (m_tileParams.imageROILeft + m_tileParams.imageROIRight);
    int ROIImgHeight = imgHeight - (m_tileParams.imageROITop + m_tileParams.imageROIBottom); 
    cout<<"Tiling.cc:: processQuadTree(): ROIImgWidth="<<ROIImgWidth<<", ROIImgHeight="<<ROIImgHeight<<endl;

    int paddedTileWidth  = m_tileParams.tileWidth + m_tileParams.paddingLeft + m_tileParams.paddingRight;
    int paddedTileHeight = m_tileParams.tileHeight + m_tileParams.paddingTop + m_tileParams.paddingBottom;
    cout<<"Tiling.cc :: processQuadTree(): paddedTileWidth="<<paddedTileWidth<<", paddedTileHeight="<<paddedTileHeight<<endl;
    float heightRatio= ROIImgHeight/(float)paddedTileHeight;
    float widthRatio = ROIImgWidth/(float)paddedTileHeight;
    float sizeRatio = heightRatio;
    if (widthRatio > sizeRatio){
      sizeRatio = widthRatio;
    }  
    cout<<"Tiling.cc:: processQuadTree(): sizeRatio="<<sizeRatio<<endl;
    
    int numPyrLevels = ceil(log2(sizeRatio));
    if (numPyrLevels < 0){
      numPyrLevels = 0;
    }
    cout<<"Tiling.cc:: processQuadTree(): numPyrLevels="<<numPyrLevels<<endl;

    int padROIImgWidth, padROIImgHeight;
    m_offsetPix.resize(2);

    
    if (pyramidMode==0){//subpyramid mode
      padROIImgWidth = ROIImgWidth;
      padROIImgHeight = ROIImgHeight;
      m_offsetPix[0] = m_tileParams.paddingLeft*pow(2.0, (double)numPyrLevels);
      m_offsetPix[1] = m_tileParams.paddingTop*pow(2.0, (double)numPyrLevels);  
    }
    
    //generates a set of multiscale tiles for an image of size imgWidth and imgHeight
    //set of non overlapping tiles, with padding, and pyramid, no offset
    //tile size is constant accross pyramid levels, but the number of tiles changes with
    //each pyramid level. at higher reolution we have more tiles than at lower res.
    //works with assembler(MSL into HiRISE) and pyramid
    if (pyramidMode==1){//pyramid mode
      padROIImgWidth = (m_tileParams.tileWidth+m_tileParams.paddingLeft+m_tileParams.paddingRight)*pow(2, (float)numPyrLevels);    
      padROIImgHeight = (m_tileParams.tileHeight+m_tileParams.paddingTop+m_tileParams.paddingBottom)*pow(2, (float)numPyrLevels);
      
      //offset from the origin of the reference
      m_offsetPix[0] = (-padROIImgWidth+ROIImgWidth)/2.0 + bBoxPixel.xl;
      m_offsetPix[1] = (-padROIImgHeight+ROIImgHeight)/2.0 + bBoxPixel.yt;
      //make sure that the offset is an integer and power of two (resampleFactor)
      m_offsetPix[0]=floor(m_offsetPix[0]/pow(2.0, numPyrLevels))*pow(2.0, numPyrLevels);
      m_offsetPix[1]=floor(m_offsetPix[1]/pow(2.0, numPyrLevels))*pow(2.0, numPyrLevels);
    }

    cout<<"Tiling.cc:: processQuadTree(): padded image="<<padROIImgWidth<<", "<<padROIImgHeight<<endl;
    cout<<"Tiling.cc:: processQuadTree(): offsetPix[0]="<<m_offsetPix[0]<<", offsetPox[1]="<<m_offsetPix[1]<<endl;

    vector<float> leftTopPixelPaddingImg;
    leftTopPixelPaddingImg.resize(2);
    
    float resampleFactor = 1;
   
    m_pyrTileParamsArray.resize(numPyrLevels+1);

    //p = 0 is the highest resolution of the reference, p=numPyrLevels is the lowest resolution, or the top of the pyramid
    for (int p = 0; p <= numPyrLevels; p++){
 
      //num tiles in resampled image
      int numHorTiles = (int)round((padROIImgWidth*resampleFactor)/paddedTileWidth);
      int numVerTiles = (int)round((padROIImgHeight*resampleFactor)/paddedTileHeight);

      leftTopPixelPaddingImg[0] = m_offsetPix[0]*resampleFactor;
      leftTopPixelPaddingImg[1] = m_offsetPix[1]*resampleFactor;
      
      cout<<"Tiling.cc::processQuadTree():: pyrLevel="<<p
	  <<", numHorTiles="<<numHorTiles<<", numVerTiles="<<numVerTiles
	  <<", leftTopPadding="<<leftTopPixelPaddingImg[0]<<", "<<leftTopPixelPaddingImg[1]
	  <<", resampleFactor="<<resampleFactor<<", padROIImgWidth="<< padROIImgWidth*resampleFactor<<endl;
	  
      //scanning the resampled image
      for (int k = 0; k <numHorTiles; k++){
	for (int l = 0; l <numVerTiles; l++){
	  
	  struct TilingParams thisTileParams;
	  thisTileParams.back_xl = k*m_tileParams.tileWidth  + leftTopPixelPaddingImg[0] - m_tileParams.paddingLeft;
	  thisTileParams.back_yt = l*m_tileParams.tileHeight + leftTopPixelPaddingImg[1] - m_tileParams.paddingTop;
	  thisTileParams.back_xr = thisTileParams.back_xl + paddedTileWidth;
	  thisTileParams.back_yb = thisTileParams.back_yt + paddedTileHeight;
	  thisTileParams.horTileIndex = k;
	  thisTileParams.verTileIndex = l;
  	  
	  m_pyrTileParamsArray[p].push_back(thisTileParams);
	  
	}
      }
      
      resampleFactor = resampleFactor*m_tileParams.pyrResampleFactor;
    }
    
  }
  /*
  void Tiling::processQuadTree(struct bbox bBoxPixel, int pyramidMode)
  {
     if (pyramidMode == 1){//pyramid mode
      cout<<"PYRAMID MODE"<<endl;
      processPyramid(bBoxPixel);
    }
    
    if (pyramidMode == 0){//subpyramid mode  
      cout<<"SUBPYRAMID MODE"<<endl;
      processSubPyramid(bBoxPixel);     
    }
  }
  */
  struct TilingParams Tiling::getBBox(int l, int hor, int ver)
  {
    struct TilingParams thisTileParams;
    
    float resampleFactor = 1.0;
    for (int i = 0; i < l; i++){
      resampleFactor = resampleFactor*m_tileParams.pyrResampleFactor;
    }

    int paddedTileWidth = m_tileParams.tileWidth+m_tileParams.paddingLeft+m_tileParams.paddingRight;
    int paddedTileHeight = m_tileParams.tileHeight+m_tileParams.paddingTop+m_tileParams.paddingBottom;
    
    vector<float> leftTopPixelPaddingImg;
    leftTopPixelPaddingImg.resize(2);
    leftTopPixelPaddingImg[0] = m_offsetPix[0]*resampleFactor;
    leftTopPixelPaddingImg[1] = m_offsetPix[1]*resampleFactor;

    thisTileParams.back_xl = leftTopPixelPaddingImg[0] + hor*m_tileParams.tileWidth - m_tileParams.paddingLeft;
    thisTileParams.back_yt = leftTopPixelPaddingImg[1] + ver*m_tileParams.tileHeight - m_tileParams.paddingTop; 
    thisTileParams.back_xr = thisTileParams.back_xl + paddedTileWidth;
    thisTileParams.back_yb = thisTileParams.back_yt + paddedTileHeight;
    thisTileParams.horTileIndex = hor;
    thisTileParams.verTileIndex = ver;
    
    leftTopPixelPaddingImg.resize(0);
    return  thisTileParams;
  }

  bool isRectangle(struct bbox bBox)
  {
    bool valid = 1;
    //cout<<"Tiling.cpp:: isRectangle(): bbox: xl="<<bBox.xl
    //  	  <<", yt="<<bBox.yt<<", xr="<<bBox.xr<<", yb="<<bBox.yb<<endl;
    if ((bBox.xl>=bBox.xr)||(bBox.yt>=bBox.yb)){
      //cout<<"Tiling.cpp:: isRectangle(): invalid bbox: xl="<<bBox.xl
      //	  <<", yt="<<bBox.yt<<", xr="<<bBox.xr<<", yb="<<bBox.yb<<endl;
      valid = 0;
    }
    
    return valid;
  }

  struct bbox computeBoxIntersection(struct bbox box1, struct bbox box2)
  {
    
    if (isRectangle(box1)==0){
      cout<<"Tiling.cpp:: computeBoxIntersection(): box1 is invalid "
	  <<"xl="<<box1.xl<<", yt="<<box1.yt<<", xr="<<box1.xr<<", yb="<<box1.yb<<endl;
    }
   
    if (isRectangle(box2)==0){
      cout<<"Tiling.cpp:: computeBoxIntersection(): box2 is invalid "
	  <<"xl="<<box2.xl<<", yt="<<box2.yt<<", xr="<<box2.xr<<", yb="<<box2.yb<<endl;
    }
   
    struct bbox intersectBox;
    
    intersectBox.xl=box1.xl;
    if (box2.xl>box1.xl){intersectBox.xl=box2.xl;}
    
    intersectBox.yt=box1.yt;
    if (box2.yt>box1.yt){intersectBox.yt=box2.yt;}
    
    intersectBox.xr=box1.xr;
    if (box2.xr<box1.xr){intersectBox.xr=box2.xr;}
    
    intersectBox.yb=box1.yb;
    if (box2.yb<box1.yb){intersectBox.yb=box2.yb;}

    return intersectBox;
  }

  void printBBox3(struct BoundingBox bbox3, string bBoxID)
  {
    cout<<bBoxID<<endl;
    cout<<"x "<<bbox3.first[0]<<", "<<bbox3.last[0]<<endl;
    cout<<"y "<<bbox3.first[1]<<", "<<bbox3.last[1]<<endl;
    cout<<"z "<<bbox3.first[2]<<", "<<bbox3.last[2]<<endl;
  }

}
