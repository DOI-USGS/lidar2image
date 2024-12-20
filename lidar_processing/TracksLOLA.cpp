/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <iomanip>
#include <math.h>

#include "../common/FileListUtils.h"
#include "../common/StringUtils.h"
#include "../common/Tiling.h"
#include "../common/ReflectanceModels.h"
#include "../geotif_processing/CoordTransform.h"
#include "TracksLOLA.h"

namespace at{

//Constructor for pointCloud which is really just a Vector3 
// with some extra data fields.
  LOLAPoint::LOLAPoint ()
  {                         
    year = 0;
    month = 0;
    day = 0;
    hour = 0;
    min = 0;
    sec = 0;
    s = 0;
  }
  
  LOLAPoint::LOLAPoint( const Eigen::Vector3f &coords,
			const int      y,
			const int      mo,
			const int      d,
			const int      h,
			const int      mi,
			const float    se,
			const int      detector)
  {                         
    year = y;
    month = mo;
    day = d;
    hour = h;
    min = mi;
    sec = se;
    s = detector;
  }
  
  bool isTimeDiff( const LOLAPoint& p, const LOLAPoint& c, const float timeThresh)
  {

    if ((p.year!=c.year) || (p.month!=c.month) || (p.day!=c.day)){
      return true;
    }
    else{
      //cout<<"TracksLOLA.cc::isTimeDiff()"<<c.sec<<", "<<p.sec<<endl;
      double duration = fabs((p.hour*3600+p.min*60+p.sec)-
			     (c.hour*3600+c.min*60+c.sec));
      //cout<<"duration = "<<duration<<endl; 
      
      if (duration <= timeThresh){
	return false;
      }
      else{
	return true;
      }
    }
  }

  
  void
  LOLAShot::init ( vector<LOLAPoint> lp,
		   vector<imgPoint>  ip,
		   vector<DEMPoint>  dp,
		   int               v,
		   float             r,
		   
		   int               ca,
		   float             fr,
		   float             fpr,
		   float             wr,
		   float             fpL,
		   float             fL)
  {
    LOLAPts = lp;
    imgPt   = ip;
    DEMPt   = dp;
    valid   = v;
    reflectance = r;
    
    //not curently used
    calc_acp = ca;
    filter_response = fr;
    featurePtRefl = fpr;
    weightRefl = wr;
    featurePtLOLA = fpL;
    filresLOLA = fL;
  }
  
  LOLAShot::LOLAShot( const vector<LOLAPoint>& pcv )
  {
    LOLAShot::init( pcv );
  }
  
  LOLAShot::LOLAShot( const LOLAPoint& pt )
  {
    vector<LOLAPoint> pcv( 1, pt );
    LOLAShot::init( pcv );
  }
  
  LOLAPoint GetPointFromIndex( const LOLAShot& shot, const int index )
  {
    const vector<LOLAPoint>& LOLAPts = shot.LOLAPts;
    
    int found = 0; 
    for( unsigned int i = 0; i < LOLAPts.size(); ++i ) {
      if( LOLAPts[i].s == index ){
	found = 1;
	return LOLAPts[i];
      }
    }
    
    if (found ==0){
      LOLAPoint noData;
      noData.x=0;
      noData.y=0;
      noData.z=0;
      
      return noData;
    }
    
  }

  int GetPointIndexFromDetector(const LOLAShot& shot, const int detectorIndex )
  {
    const vector<LOLAPoint>& LOLAPts = shot.LOLAPts;
      
    int found = 0;
    for( unsigned int i = 0; i < LOLAPts.size(); ++i ) {
      if( LOLAPts[i].s == detectorIndex ){
	found = 1;
	return i;
      }
    }
    if (found == 0){
      return (-1);
    }
  }
    
  
  vector<vector<LOLAShot> > LOLAFileRead( const string& f )
  {
    
    // This is specifically for reading the RDR_*PointPerRow_csv_table.csv files only.
    // It will ignore lines which do not start with year (or a value that can be converted
    // into an integer greater than zero, specfically).
    
    char temp[255];
    ifstream file( f.c_str() );
    if( !file ) {
      cout<<"Unable to open track file "<<f<<endl;
    }
 

    vector<vector<LOLAShot> > LOLATracks(1);
    
    // Establish 3d bounding box for acceptable LOLA data. 
    BoundingBox LOLArange;
    LOLArange.first.resize(3);LOLArange.last.resize(3);
    LOLArange.first[0]= -180; LOLArange.last[0]=  360;
    LOLArange.first[1]=  -90; LOLArange.last[1]=   90;
    LOLArange.first[2]= 1720; LOLArange.last[2]= 1750;
    
    int trackIndex = 0;
    LOLAPoint currPt = LOLAPoint();
    LOLAPoint prevPt = LOLAPoint();
    
    string line;
    while( getline(file, line, '\n') ) {
      
      strncpy(temp, line.c_str(), 255);
      const char* token = strtok(temp, ",");
      
      int ret = sscanf(token, "%d-%d-%dT%d:%d:%g", &currPt.year, &currPt.month, &currPt.day, &currPt.hour, &currPt.min, &currPt.sec);
      if( currPt.year <= 0 ) { continue; }
      
      token = strtok(NULL, ",");
      ret += sscanf(token, "%lg", &currPt.x);
      
      token = strtok(NULL, ",");
      ret = ret + sscanf(token, "%lg", &currPt.y);
      
      token = strtok(NULL, ",");
      ret = ret + sscanf(token, "%lg", &currPt.z);
  
      for (int i = 0; i < 8; i++)
	token = strtok(NULL, ",");
      ret = ret + sscanf(token, "%d", &currPt.s);
      
      if (ret != 10){
	fprintf(stderr, "Failed to read line %s.\n", line.c_str());
	cout<<"Failed to read line"<<endl;
      }
      
      if ((currPt.x > LOLArange.first[0]) && (currPt.x < LOLArange.last[0]) &&
	  (currPt.y > LOLArange.first[1]) && (currPt.y < LOLArange.last[1]) &&
	  (currPt.z > LOLArange.first[2]) && (currPt.z < LOLArange.last[2])){
	
	if (LOLATracks[trackIndex].empty()) {
	    LOLATracks[trackIndex].push_back( LOLAShot(currPt) );
	}
	else {
	  if( isTimeDiff(prevPt, currPt, 3000) ) { //new track
	    trackIndex++;
	    LOLATracks.resize(trackIndex+1); //start new track
	    LOLATracks[trackIndex].push_back( LOLAShot(currPt) ); //put the new shot in it
	  }
	  else { //same track
	    if(isTimeDiff(prevPt, currPt, 0) ) { //new shot
	       LOLATracks[trackIndex].push_back( LOLAShot(currPt) );
	    }
	    else { //same shot
	      LOLATracks[trackIndex].back().LOLAPts.push_back( currPt );
	    }
	  }
	}
	
	//copy current pc into prevPt
	prevPt = currPt;
      }
      else{
	cout<<"TracksLOLA.cc::LOLAFileRead "<<currPt.x<<", "<<currPt.y<<", "<<currPt.z<<endl;
      }
    } 
    file.close();

  return LOLATracks; 

  }

  vector<vector<LOLAShot> > LOLAFilesRead( const string& LOLAFilename )
  {
    std::vector<std::string> trackFilenames = AccessDataFilesFromInput(LOLAFilename,"csv");
    cout<<"num trackFilenames = "<<trackFilenames.size()<<endl;
    vector<vector<LOLAShot> > LOLATracks;
    for (unsigned int i = 0; i < trackFilenames.size(); i++){
      cout<<trackFilenames[i]<<endl; 
      vector<vector<LOLAShot> > LOLATracksInFile = LOLAFileRead(trackFilenames[i]);
      LOLATracks.reserve(LOLATracks.size() + distance(LOLATracksInFile.begin(),LOLATracksInFile.end()));
      LOLATracks.insert(LOLATracks.end(), LOLATracksInFile.begin(), LOLATracksInFile.end());
    }
    return LOLATracks;
  }
  
  int ComputeAllReflectance( vector< vector<LOLAShot> >& LOLATracks,  
			     const Eigen::Vector3f& cameraPosition, 
			     const Eigen::Vector3f& lightPosition)
  {
    int numValidReflPts = 0;
    int numTotalPts = 0;
    
    for( unsigned int k = 0; k < LOLATracks.size(); ++k ){
      for( unsigned int i = 0; i < LOLATracks[k].size(); ++i ){
	
	numTotalPts++;
	
	LOLATracks[k][i].valid = 0;
	LOLATracks[k][i].reflectance = -1;
	
	if ((LOLATracks[k][i].LOLAPts.size()>0) && (LOLATracks[k][i].LOLAPts.size()<=5)){
	  
	  //select center point
	  //LOLAPoint centerPt = GetPointFromIndex( LOLATracks[k][i].LOLAPts, 1);
	  LOLAPoint centerPt = GetPointFromIndex( LOLATracks[k][i], 1);
	  if ((centerPt.x!=0)||(centerPt.y!=0)||(centerPt.z!=0)){
	    Eigen::Vector3f lonLatRad;
	    lonLatRad[0] = M_PI*centerPt.x/180.0;
	    lonLatRad[1] = M_PI*centerPt.y/180.0;
	    lonLatRad[2] = 1000*centerPt.z;
	    
	    Eigen::Vector3f center;
	    center = geographicToCartesian(lonLatRad);
	    
	    Eigen::Vector3f points[4];
	    bool validPoints[4] = {true, true, true, true};
	    
	    // compute average of normal for each triangle,
	    // LOLA shots are positioned like number 5 of a die
	    // 1 is center dot, 2 is left, 3 is top, 4 is right, 5 is bottom
	    for (int m = 2; m <= 5; m++){
	      
	      //LOLAPoint pt = GetPointFromIndex( LOLATracks[k][i].LOLAPts, m);
	      LOLAPoint pt = GetPointFromIndex( LOLATracks[k][i], m);
	      if ((pt.x!=0)||(pt.y!=0)||(pt.z!=0)){
		Eigen::Vector3f lonLatRad;
		lonLatRad[0] = M_PI*pt.x/180.0;
		lonLatRad[1] = M_PI*pt.y/180.0;
		lonLatRad[2] = 1000*pt.z;
		points[m-2]  = geographicToCartesian(lonLatRad);
	      }
	      else{
		validPoints[m-2] = false;
	      }
	    }
	    
	    Eigen::Vector3f normal(0.0, 0.0, 0.0);
	    
	    for (int m = 0; m < 4; m++){
	      int n;
	      if (m == 3){n = 0;}
	      else {n = m + 1;}	
	      if (validPoints[m] && validPoints[n]){
		normal += ComputeNormalFrom3DPointsGeneral(center, points[m], points[n]);
	      }
	    }
	    
	    if (normal[0]!=0.0 || normal[1]!=0.0 || normal[2]!=0.0){
	      
	      normal = normal.normalized();	    
	      
	      center[0] = center[0]/1000.0;
	      center[1] = center[1]/1000.0;
	      center[2] = center[2]/1000.0;
	      
	      //cout<<"TracksLOLA.cc::ComputeAllReflectance(): center="
	      //<<center[0]<<", "<<center[1]<<", "<<center[2]<<endl;
	      //cout<<"TracksLOLA.cc::ComputeAllReflectance(): camPos="
	      //<<cameraPosition[0]<<", "<<cameraPosition[1]<<", "<<cameraPosition[2]<<endl;
	      //cout<<"TracksLOLA.cc::ComputeAllReflectance(): sunPos="
	      //<<lightPosition[0]<<", "<<lightPosition[1]<<", "<<lightPosition[2]<<endl;
	      //cout<<"TracksLOLA.cc::ComputeAllReflectance(): normal="
	      //<<normal[0]<<", "<<normal[1]<<", "<<normal[2]<<endl;
	      //cout<<"----------------------------------------------------"<<endl;
	      LOLATracks[k][i].reflectance = ComputeLunarLambertianReflectanceFromNormal(lightPosition, cameraPosition, center, normal);	      
	      numValidReflPts++;
	      //shots for which a valid reflectance can be computed
	      LOLATracks[k][i].valid = 1;
	    }
	  }
	}
      }//i
    }//k
    
    cout<<"TracksLOLA.cc:: ComputeAllReflectance(): numValidReflPts = "<<numValidReflPts<<", numTotalPts="<<numTotalPts<<endl;
    
    return numValidReflPts;
  }
  

  void SaveReflectancePoints( const vector<vector<LOLAShot> >& LOLATracks, 
			      const Eigen::Vector2f&           gain_bias,
			      const string&                    filename)
  {
    //boost::filesystem::path p( filename );
    for( unsigned int k = 0; k < LOLATracks.size(); ++k ){
      //ostringstream os;
      //os << p.stem().string() << "_" << k << ".txt";

      string outFilename = GetFilenameNoPath(GetFilenameNoExt(filename))+"_"+NumberToString(k)+".txt";
      
      //ofstream file( os.str().c_str() );
      ofstream file(outFilename.c_str());
      if( !file ) {
	//cout<< "Can't open reflectance output file " << os.str()<<endl;
	cout<< "Can't open reflectance output file " << outFilename<<endl;
      }
      
      for( unsigned int i = 0; i < LOLATracks[k].size(); ++i){
	if ( (LOLATracks[k][i].valid == 1) && 
	     (LOLATracks[k][i].reflectance > 0) ){//valid track and non-zero reflectance
	  file << ( gain_bias(0)*LOLATracks[k][i].reflectance + gain_bias(1) ) << endl;
	}
	else{ file << "-1" << endl; }
      }
      file.close();
    }
  }
  
  //saves the image points corresponding to a detectNum
  void SaveImagePoints( const vector<vector<LOLAShot> >& LOLATracks,
			const int                        detectNum,
			const string&                    filename)
  {
    //boost::filesystem::path p( filename );
    for (unsigned int k = 0; k < LOLATracks.size(); ++k ){
      //ostringstream os;
      //os << p.stem().string() << "_" << k << ".txt";
      //ofstream file( os.str().c_str() );
      
      string outFilename = GetFilenameNoPath(GetFilenameNoExt(filename)) + "_" + NumberToString(k) + ".txt";
      ofstream file(outFilename.c_str() );
 
      if( !file ) {
	//cout<< "Can't open image point output file " << os.str()<<endl;
	cout<< "Can't open reflectance output file " << outFilename<<endl;
      }
      
      for( unsigned int i = 0; i < LOLATracks[k].size(); ++i){
	bool found = false;
	for( unsigned int j = 0; j < LOLATracks[k][i].LOLAPts.size(); ++j){
	  if( (LOLATracks[k][i].LOLAPts[j].s == detectNum) && 
	      (LOLATracks[k][i].valid       == 1)           ){
	    found = true;
	    file <<  LOLATracks[k][i].imgPt[j].val << endl;
	  }
	}
	if( !found ){ file << "-1" << endl; }
      }
      file.close();
    }
  }
  
  //saves the altitude  points (LOLA) corresponding to a sensor ID = detectNum
  void SaveAltitudePoints( const vector<vector<LOLAShot> >& LOLATracks,
			   const int                        detectNum,
			   const string&                    filename)
  {
    //boost::filesystem::path p( filename );
    for (unsigned int t = 0; t < LOLATracks.size(); ++t ){
      //ostringstream os;
      //os << p.stem().string() << "_" << t << ".txt";
      //ofstream file( os.str().c_str() );
      
      string outFilename = GetFilenameNoPath(GetFilenameNoExt(filename)) + "_" + NumberToString(t) + ".txt";
      ofstream file(outFilename.c_str() );
   
      if( !file ) {
	//cout<< "Can't open altitude output file " << os.str()<<endl;
	cout<< "Can't open reflectance output file " << outFilename<<endl;
      }
      
      for( unsigned int s = 0; s < LOLATracks[t].size(); ++s ){
	bool found = false;
	for( unsigned int j = 0; j < LOLATracks[t][s].LOLAPts.size(); ++j ){
	  if( (LOLATracks[t][s].LOLAPts[j].s == detectNum) && 
	      (LOLATracks[t][s].valid       == 1)           ){
	    found = true;
	    file << LOLATracks[t][s].LOLAPts[j].z << endl;
	  }
	}
	if( !found ){ file << "-1" << endl; }
      }
      file.close();
    }
  }


  
  //computes the min and max lon and lat of the LOLA data
  BoundingBox FindShotBounds(const vector<vector<LOLAShot> >& LOLATracks) 
  {
    Eigen::Vector2f min( std::numeric_limits<float>::max(), 
			 std::numeric_limits<float>::max() );
    Eigen::Vector2f max( std::numeric_limits<float>::min(), 
			 std::numeric_limits<float>::min() );
    for(     unsigned int i = 0; i < LOLATracks.size();              ++i){
      for(   unsigned int j = 0; j < LOLATracks[i].size();           ++j){
	for( unsigned int k = 0; k < LOLATracks[i][j].LOLAPts.size(); ++k){
	  double lon = LOLATracks[i][j].LOLAPts[k].x;
	  double lat = LOLATracks[i][j].LOLAPts[k].y;
	  
	  if (lat < min.y() ){ min.y() = lat; }
	  if (lon < min.x() ){ min.x() = lon; }
	  if (lat > max.y() ){ max.y() = lat; }
	  if (lon > max.x() ){ max.x() = lon; }
	}
      }
    }
    
    BoundingBox bounds;
    bounds.first.resize(2);
    bounds.last.resize(2);
    bounds.first[0]=min.x();  bounds.last[0]=max.x();
    bounds.first[1]=min.y();  bounds.last[1]=max.y();
    
    return bounds;
  }
  
  Eigen::Vector4f FindMinMaxLat( const vector<vector<LOLAShot> >& LOLATracks )
  {
    BoundingBox bounds = FindShotBounds( LOLATracks );
    Eigen::Vector4f coords;
    coords(0) = bounds.first[1];
    coords(1) = bounds.last[1];
    coords(2) = bounds.first[0];
    coords(3) = bounds.last[0];
    
    return coords;
  }
    
   void DrawTracksOnImage(vector<vector<LOLAShot> > tracks,
			  GDALDataset* imageDS, string trackType, int thikness)
  {

    int imageWidth  = GDALGetRasterXSize(imageDS);
    int imageHeight = GDALGetRasterYSize(imageDS);
    GDALRasterBand* imageBand = imageDS->GetRasterBand(1);
    float noDataValue = imageBand->GetNoDataValue();
    
    float *image = new float[imageWidth*imageHeight];
    
    imageBand->RasterIO(GF_Read, 0, 0, imageWidth, imageHeight,  
			image, imageWidth, imageHeight, 
			GDT_Float32, 0, 0);

    float minImageVal = 1.0;//imageBand->GetMinimum();
    float maxImageVal = 0.0;//imageBand->GetMaximum();
    for (int row = 0; row < imageHeight; row++){
      for (int col = 0; col < imageWidth; col++){
	if (image[row*imageWidth+col]!=noDataValue){
	  if (image[row*imageWidth+col] > maxImageVal){
	    maxImageVal = image[row*imageWidth+col];
	  }
	  if (image[row*imageWidth+col] < minImageVal){
	    minImageVal = image[row*imageWidth+col];
	  }
        }
      }
    }
    

    double gainDisplay = 1.0/(maxImageVal-minImageVal);
    double biasDisplay = -gainDisplay*minImageVal;
    
    cout<<"minImageVal="<<minImageVal<<", maxImageVal="<<maxImageVal<<", gain="<<gainDisplay<<", bias="<<biasDisplay<<endl;
    
    for (int row = 0; row < imageHeight; row++){
      for (int col = 0; col < imageWidth; col++){
	image[row*imageWidth+col] = gainDisplay*image[row*imageWidth+col]+biasDisplay;
      }
    }
    
    //draw the aligned tracks point
    for (unsigned int t = 0; t < tracks.size(); t++){ //for each track
      for (unsigned int s = 0; s < tracks[t].size(); s++){//for each shot	 
        if ((tracks[t][s].valid == 1)){

	  int index = GetPointIndexFromDetector(tracks[t][s], 1);
	  
	  //display the center point if it is valid
	  if (index!=-1){
	    int ptRow = (int)floor(tracks[t][s].imgPt[index].y);
	    int ptCol = (int)floor(tracks[t][s].imgPt[index].x);
	    
	    int xl = ptCol - thikness; if(xl < 0) {xl = 0;}
	    int xr = ptCol + thikness; if(xr > imageWidth-1) {xr = imageWidth-1;}
	    int yt = ptRow - thikness; if(yt < 0) {yt = 0;}
	    int yb = ptRow + thikness; if(yb > imageHeight-1) {yb = imageHeight-1;}
	    
	    if (trackType.compare("tracks")==0){
	      for (int row = yt; row < yb; row++){
		for (int col = xl; col < xr; col++){
		  image[row*imageWidth+col] = 1.0;
		}
	      }
	    }
	    if (trackType.compare("reflectance")==0){
	      for (int row = yt; row < yb; row++){
		for (int col = xl; col < xr; col++){
		  if ((col==xl) || (col==xr-1)){
		    image[row*imageWidth+col] = 1.0;
		  }
		  else{
		    image[row*imageWidth+col] = gainDisplay*tracks[t][s].reflectance+biasDisplay;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    
    imageBand->RasterIO(GF_Write, 0, 0, imageWidth, imageHeight,  
			image, imageWidth, imageHeight, 
			GDT_Float32, 0, 0);    
  }
  
  vector<string> trackSplitter(string lidarFilename, string resultsDirname)
  {
    vector<string> outFilenames;
    
    ifstream file(lidarFilename.c_str() );
    if( !file ) {
      cerr << "Unable to open lidar file " << lidarFilename << endl;
      return outFilenames;
    }
    
    LOLAPoint currPt;
    LOLAPoint prevPt;
    string line;
    ofstream outfile;
    
    while( getline(file, line, '\n') ) {
      
      vector<string> elements = split((const string)line, ',');
      vector<string> timeAndDate = split (elements[0], 'T');
      vector<string> ymd = split(timeAndDate[0], '-');
      vector<string> hms = split(timeAndDate[1], ':');
      
      if ((ymd.size()==3) && (hms.size()==3)){
	currPt.year = atoi(ymd[0].c_str());
	currPt.month = atoi(ymd[1].c_str());
	currPt.day = atoi(ymd[2].c_str());
	currPt.hour = atoi(hms[0].c_str());
	currPt.min = atoi(hms[1].c_str());
	currPt.sec = atof(hms[2].c_str());
      }
      
      if( currPt.year <= 0 ) { continue; }
      
      if( prevPt.year == 0 || isTimeDiff(prevPt, currPt, 3000)) { //new track
	
	if( outfile.is_open() ){
	  outfile.close();
	}
	
	string outname=resultsDirname + "/" + GetFilenameNoPath(GetFilenameNoExt(lidarFilename))
	  + "_" + NumberToString(currPt.year)+"_"+
	  NumberToString(currPt.month)+ "_" +
	  NumberToString(currPt.day)+".csv";
	
	outFilenames.push_back(outname);   
	
	//open outputfile if resultsDirname exists
	if (resultsDirname.compare("")!=0){
	  outfile.open(outname.c_str());
	}
	
      }
      //write outputfile if resultsDirname exists
      if (resultsDirname.compare("")!=0){
	outfile << line << endl;
      }
      
      prevPt.year=currPt.year;
      prevPt.month=currPt.month;
      prevPt.day=currPt.day;
      prevPt.hour=currPt.hour;
      prevPt.min=currPt.min;
      prevPt.sec=currPt.sec;
      
    } 
    
    outfile.close();
    file.close();
    
    for (int i = 0; i < outFilenames.size(); i++){
      cout<<outFilenames[i]<<endl;
    }
    
    return outFilenames;
  }

  
  int test(vector<string>outTrackFilenames, string testTrackListFilename)
  {
    int testFailed = 0;
    
    vector<string> testTrackFilenames;
    
    ifstream file(testTrackListFilename.c_str() );
    if( !file ) {
      cerr << "Unable to open testTrackListFilename file " << testTrackListFilename << endl;
      testFailed = 1;
      
    }
    else{  
      string line;
      while( getline(file, line, '\n') ) {
	testTrackFilenames.push_back(line); 
      }
      
      if (testTrackFilenames.size()!=outTrackFilenames.size()){
	cout<<"incorrect number of files:"<< outTrackFilenames.size()<<" instead of "<<testTrackFilenames.size()<<endl;
	testFailed = 1;
      }
      else{
	for (int i = 0; i < testTrackFilenames.size(); i++){
	  if (testTrackFilenames[i].compare(outTrackFilenames[i])!=0){
	    testFailed = 1;
	    cout<<"incorrect line "<<outTrackFilenames[i]<<endl;
	  }
	}
	
      }
      
    }
    
    return testFailed;
  }
}
