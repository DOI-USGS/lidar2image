#ifndef ISIS_INTERFACE_H
#define ISIS_INTERFACE_H

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

//#if 0
#include "SurfacePoint.h"
#include "Latitude.h"
#include "Longitude.h"
//#endif  

using namespace std;
using namespace Isis;
namespace at
{

// A program demonstrating how to call point_to_pixel(), pixel_to_vector(),
// and camera_center() for an ISIS frame camera.

// Very important note: This code assumes that pixel indices start from 0,
// not from 1 like in ISIS.

// Utilities
string print_vec(vector<double> const& vec);
vector<double> diff(vector<double> const& vec1, vector<double> const& vec2);
vector<double> prod(double val, vector<double> const& vec);
vector<double> div(vector<double> const& vec, double val);
vector<double> normalize(vector<double> const& vec);

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
  virtual vector<double>
  camera_rotation() const = 0;
  Isis::Camera* GetCamera(){return m_camera;}
  
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
  virtual vector<double>
  camera_rotation() const;
protected:
  
  // Custom Variables
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  Isis::CameraDetectorMap   *m_detectmap;
  mutable Isis::AlphaCube   m_alphacube;
  
  vector<double> m_center;
  vector<double> m_rotation;
};

//#if 0

// ISIS linescan camera class
class IsisInterfaceLineScan : public IsisInterface {

public:
  IsisInterfaceLineScan( std::string const& filename );

  virtual ~IsisInterfaceLineScan() {}

  virtual std::string type()  { return "LineScan"; }

  // Standard Methods
  //-------------------------------------------------
  
  // Note: The smallest pixel is (0, 0), not (1, 1)!
  virtual vector<double>
  point_to_pixel( vector<double> const& point ) const;
  virtual vector<double>
  pixel_to_vector( vector<double> const& pix ) const;
  virtual vector<double>
  camera_center() const;
  virtual vector<double>
  camera_rotation() const;
protected:

  // Custom Variables
  Isis::CameraDistortionMap *m_distortmap;
  Isis::CameraFocalPlaneMap *m_focalmap;
  Isis::CameraDetectorMap   *m_detectmap;
  mutable Isis::AlphaCube   m_alphacube; // Doesn't use const
};
  
//#endif
  
  }
#endif /* ISIS_INTERFACE_H */
