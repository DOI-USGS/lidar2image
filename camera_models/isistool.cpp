//#include "Isis.h"
#include "Camera.h"
#include "CameraFactory.h"
#include "TProjection.h"
#include "ProjectionFactory.h"
#include "ProcessRubberSheet.h"
#include "IException.h"
#include "Pvl.h"
#include "IString.h"
#include "PushFrameCameraDetectorMap.h"
#include "Target.h"
#include "CameraDistortionMap.h"
#include "CameraFocalPlaneMap.h"
#include <vector>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;
using namespace Isis;

// A program demonstrating how to call point_to_pixel(), pixel_to_vector(),
// and camera_center() for an ISIS frame camera.

// Very important note: This code assumes that pixel indices start from 0,
// not from 1 like in ISIS.

// Utilities
string print_vec(vector<double> const& vec){
  ostringstream oss;
  oss.precision(18);
  oss << "(";
  for (int s = 0; s < (int)vec.size()-1; s++){
    oss << vec[s] << ", ";
  }
  if (!vec.empty()) oss << vec[vec.size()-1];
  oss << ")";
  return oss.str();
}

vector<double> diff(vector<double> const& vec1, vector<double> const& vec2){
  vector<double> ans = vec1;
  for (int i = 0; i < (int)vec1.size(); i++){
    ans[i] -= vec2[i];
  }
  return ans;
}

vector<double> prod(double val, vector<double> const& vec){
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] *= val;
  }
  return ans;
}

vector<double> div(vector<double> const& vec, double val){
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] /= val;
  }
  return ans;
}

vector<double> normalize(vector<double> const& vec){
  double len = 0;
  for (int i = 0; i < (int)vec.size(); i++){
    len += vec[i]*vec[i];
  }
  len = sqrt(len);
  vector<double> ans = vec;
  for (int i = 0; i < (int)vec.size(); i++){
    ans[i] /= len;
  }
  return ans;
}

// Abstract base class
class IsisInterface {
public:
  IsisInterface( string const& file );
  virtual ~IsisInterface(); // Can't be declared here since we have
  // incomplete types from Isis.

  virtual string type() = 0;
  static IsisInterface* open( string const& filename );

  // Standard Methods
  //------------------------------------------------------

  // These are the standard camera request; IsisInterface allows for
  // them to be customized for the type of camera so that they are
  // fast and not too full of conditionals.

  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const = 0;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const = 0;
  virtual vector<double>
  camera_center() const = 0;

protected:
  // Standard Variables
  //------------------------------------------------------
  Isis::Pvl * m_label;
  Isis::Camera * m_camera;
  Isis::Cube * m_cube;

};

// ISIS frame camera class
class IsisInterfaceFrame : public IsisInterface {
  
public:
  IsisInterfaceFrame( string const& filename );
  
  virtual string type()  { return "Frame"; }
  
  // Standard Methods
  //-------------------------------------------------
  
  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const;
  virtual vector<double>
  camera_center() const;
  
protected:
  
  // Custom Variables
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  Isis::CameraDetectorMap   *m_detectmap;
  mutable Isis::AlphaCube   m_alphacube;
  
  vector<double> m_center;
};

IsisInterface::IsisInterface( string const& file ) {
  // Opening labels and camera
  Isis::FileName ifilename( QString::fromStdString(file) );
  m_label = new Isis::Pvl();
  m_label->read( ifilename.expanded() );

  // Opening Isis::Camera
  m_cube = new Isis::Cube(QString::fromStdString(file));
  m_camera = Isis::CameraFactory::Create( *m_cube );
}

IsisInterface::~IsisInterface(){
  delete m_label;
  delete m_camera;
  delete m_cube;
}

IsisInterface* IsisInterface::open( string const& filename ) {
  IsisInterface* result;

  // Opening Labels (This should be done somehow though labels)
  Isis::FileName ifilename( QString::fromStdString(filename) );
  Isis::Pvl label;
  label.read( ifilename.expanded() );

  Isis::Cube tempCube(QString::fromStdString(filename));
  Isis::Camera* camera = Isis::CameraFactory::Create( tempCube );

  switch ( camera->GetCameraType() ) {
  case 0:
    // Framing Camera
    if ( camera->HasProjection() ){
      throw string("Don't support Isis Camera Type with projection at this moment");
      //result = new IsisInterfaceMapFrame( filename );
    }
    else{
      result = new IsisInterfaceFrame( filename );
    }
    break;
  case 2:
    // Linescan Camera
    if ( camera->HasProjection() ){
      throw string("Don't support Projected Linescan Isis Camera Type at this moment");
      //result = new IsisInterfaceMapLineScan( filename );
    }
    else{
      throw string("Don't support Linescan Isis Camera Type at this moment");
      //result = new IsisInterfaceLineScan( filename );
    }
    break;
  default:
    throw string("Don't support unknown Isis Camera Type at this moment");
  }
  return result;
}

IsisInterfaceFrame::IsisInterfaceFrame( string const& filename ) :
  IsisInterface(filename), m_alphacube( *m_cube ) {

  // Gutting Isis::Camera
  m_distortmap = m_camera->DistortionMap();
  m_focalmap   = m_camera->FocalPlaneMap();
  m_detectmap  = m_camera->DetectorMap();

  // Calculating Center (just once)
  m_center.resize(3);
  m_camera->instrumentPosition(&m_center[0]);
  for (int i = 0; i < (int)m_center.size(); i++){
    m_center[i] *= 1000; // convert to meters
  }
}

vector<double>
IsisInterfaceFrame::point_to_pixel( vector<double> const& point ) const{

  // Note: The smallest pixel is (0, 0), not (1, 1)!

  vector<double> look = normalize( diff(point, m_center) );

  vector<double> lookB_copy(3);
  std::copy( look.begin(), look.end(), lookB_copy.begin() );
  lookB_copy = m_camera->bodyRotation()->J2000Vector(lookB_copy);
  lookB_copy = m_camera->instrumentRotation()->ReferenceVector(lookB_copy);
  std::copy( lookB_copy.begin(), lookB_copy.end(), look.begin() );
  look = div(prod(m_camera->FocalLength(), look), fabs(look[2]) );

  // Back Projecting
  m_distortmap->SetUndistortedFocalPlane( look[0], look[1] );
  m_focalmap->SetFocalPlane( m_distortmap->FocalPlaneX(),
                             m_distortmap->FocalPlaneY() );
  m_detectmap->SetDetector( m_focalmap->DetectorSample(),
                            m_focalmap->DetectorLine() );

  vector<double> ans(2);
  ans[0] = m_alphacube.BetaSample( m_detectmap->ParentSample() ) - 1;
  ans[1] = m_alphacube.BetaLine( m_detectmap->ParentLine() ) - 1;
  return ans;
}

vector<double> 
IsisInterfaceFrame::pixel_to_vector( vector<double> const& px ) const {

  // Note: The smallest pixel is (0, 0), not (1, 1)!

  m_detectmap->SetParent( m_alphacube.AlphaSample(px[0]+1),
                          m_alphacube.AlphaLine(px[1]+1));
  m_focalmap->SetDetector( m_detectmap->DetectorSample(),
                           m_detectmap->DetectorLine() );
  m_distortmap->SetFocalPlane( m_focalmap->FocalPlaneX(),
                               m_focalmap->FocalPlaneY() );
  vector<double> look(3);
  look[0] = m_distortmap->UndistortedFocalPlaneX();
  look[1] = m_distortmap->UndistortedFocalPlaneY();
  look[2] = m_distortmap->UndistortedFocalPlaneZ();

  look = normalize( look );
  look = m_camera->instrumentRotation()->J2000Vector(look);
  look = m_camera->bodyRotation()->ReferenceVector(look);
  return look;
}

vector<double> IsisInterfaceFrame::camera_center() const{
  return m_center;
}

int main( int argc, char *argv[] ) {

  for (int s = 0; s < argc; s++) cout << argv[s] << " ";
  cout << endl;

  if (argc < 2){
    std::cerr << "Usage: " << argv[0] << " cubFile.cub"  << std::endl;
  }

  // Load the cube
  string cubFile = argv[1];
  IsisInterfaceFrame cam(cubFile);

  // An xyz point
  vector<double> xyz;
  xyz.push_back(1.55403e+06);
  xyz.push_back(86956.4);
  xyz.push_back(767224);

  vector<double> pix = cam.point_to_pixel(xyz);
  vector<double> vec = cam.pixel_to_vector(pix);
  vector<double> ctr = cam.camera_center();

  std::cout << "xyz point:       " << print_vec(xyz) << std::endl;
  std::cout << "point_to_pixel:  " << print_vec(pix) << std::endl;
  std::cout << "pixel_to_vector: " << print_vec(vec) << std::endl;
  std::cout << "camera center:   " << print_vec(ctr) << std::endl;
}
