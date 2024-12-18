//Eigen includes
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "IsisInterfaceATK.h"


using namespace std;
using namespace Isis;

namespace at
{
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
	result = new IsisInterfaceLineScan( filename );
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
    //Calculating Rotation (just once)
    std::vector<double> rot_inst = m_camera->instrumentRotation()->Matrix();
    std::vector<double> rot_body = m_camera->bodyRotation()->Matrix();
    Eigen::Matrix3d RotInst;
    RotInst(0,0)=rot_inst[0];   RotInst(0,1)=rot_inst[1];   RotInst(0,2)=rot_inst[2];
    RotInst(1,0)=rot_inst[3];   RotInst(1,1)=rot_inst[4];   RotInst(1,2)=rot_inst[5];
    RotInst(2,0)=rot_inst[6];   RotInst(2,1)=rot_inst[7];   RotInst(2,2)=rot_inst[8];
    Eigen::Matrix3d RotBody;
    RotBody(0,0)=rot_body[0];   RotBody(0,1)=rot_body[1];   RotBody(0,2)=rot_body[2];
    RotBody(1,0)=rot_body[3];   RotBody(1,1)=rot_body[4];   RotBody(1,2)=rot_body[5];
    RotBody(2,0)=rot_body[6];   RotBody(2,1)=rot_body[7];   RotBody(2,2)=rot_body[8];
    Eigen::Matrix3d Rot = RotBody*(RotInst.transpose());
    m_rotation.resize(9);
    m_rotation[0]=Rot(0,0);    m_rotation[1]=Rot(0,1);   m_rotation[2]=Rot(0,2);
    m_rotation[3]=Rot(1,0);    m_rotation[4]=Rot(1,1);   m_rotation[5]=Rot(1,2);
    m_rotation[6]=Rot(2,0);    m_rotation[7]=Rot(2,1);   m_rotation[8]=Rot(2,2);
    //= m_camera->bodyRotation()->Matrix();
    //R_body*transpose(R_inst)
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
  vector<double> IsisInterfaceFrame::camera_rotation() const{
    return m_rotation;
  }

  // IsisInterfaceLineScan class implementation
  IsisInterfaceLineScan::IsisInterfaceLineScan( std::string const& filename ):
    IsisInterface(filename), m_alphacube( *m_cube ) {
    
    // Gutting Isis::Camera
    m_distortmap = m_camera->DistortionMap();
    m_focalmap   = m_camera->FocalPlaneMap();
    m_detectmap  = m_camera->DetectorMap();  
  }
  
  vector<double>
  IsisInterfaceLineScan::point_to_pixel( vector<double> const& point ) const
  {
    using namespace Isis;
    vector<double> pix(2);
    
  // Note: The smallest pixel is (0, 0), not (1, 1)!
  SurfacePoint surfP(Displacement(point[0], Displacement::Meters),
    Displacement(point[1], Displacement::Meters),
    Displacement(point[2], Displacement::Meters));
  
  if(m_camera->SetUniversalGround(surfP.GetLatitude().degrees(),
				  surfP.GetLongitude().degrees(),
				  surfP.GetLocalRadius().meters())){
    pix[0] = m_camera->Sample()-1;
    pix[1] = m_camera->Line()-1;
    //cout<<"pix="<<pix[0]<<", "<<pix[1]<<endl;
   }
  else{
    cout<<"Could not project point into the camera."<<endl;
    pix[0]=0;
    pix[1]=0;
    //cout<<"lat="<<surfP.GetLatitude().degrees()<<endl;
    //cout<<"lon="<<surfP.GetLongitude().degrees()<<endl;
    //cout<<"rad="<<surfP.GetLocalRadius().meters()<<endl;
  }
  return pix;
 }
  
vector<double> 
IsisInterfaceLineScan::pixel_to_vector( vector<double> const& px ) const {

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

vector<double> IsisInterfaceLineScan::camera_center() const{
  throw string("Computation of camera center is not implemented")
    + " for linescan cameras.";
}
vector<double> IsisInterfaceLineScan::camera_rotation() const{
  throw string("Computation of camera rotation is not implemented")
    + " for linescan cameras.";
}
}

