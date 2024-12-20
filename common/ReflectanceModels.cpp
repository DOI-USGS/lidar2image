/******************************************************************************
 * NOTE: Original code authored by Ara Nefien in 2006-2007 and placed in the 
 * public domain per U.S. government policy and personal communication.  
 ******************************************************************************
 * 
 * SPDX-License-Identifier: CC0-1.0
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "ReflectanceModels.h"
using namespace std;

namespace at{

/* Deprecated?  Not used anywhere
Vector3 ComputePlaneNormalFrom3DPoints(vector<Vector3> pointArray)
{
  Vector3 normal;
  Matrix<float,5,4> rhs;
  Vector<float,3> lhs;
  for (unsigned int i = 0; i < pointArray.size(); i++){
    rhs(i,0)=pointArray[i][0];
    rhs(i,1)=pointArray[i][1]; 
    rhs(i,2)=pointArray[i][2];
    rhs(i,3)=1;
  }
  
  solve_symmetric_nocopy(rhs,lhs);
  return normal;

}
*/

//Vector3 ComputeNormalFrom3DPointsGeneral(Vector3 p1, Vector3 p2, Vector3 p3) {
//  return -normalize(cross_prod(p2-p1,p3-p1));
//}

Eigen::Vector3f ComputeNormalFrom3DPointsGeneral(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3) 
{
  Eigen::Vector3f cp=-((p2-p1).cross(p3-p1));
  return cp.normalized();       
  //return -normalize(cross_prod(p2-p1,p3-p1));
}

float ComputeLunarLambertianReflectanceFromNormal(Eigen::Vector3f sunPos, Eigen::Vector3f viewPos, 
						  Eigen::Vector3f xyz, Eigen::Vector3f normal) 
{
  float reflectance;
  float L;

  if (isnan(normal[0])||(isnan(normal[1])) || (isnan(normal[2]))){
    return 0.0;
  }
  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  //sun coordinates relative to the xyz point on the Moon surface
  //Vector3 sunDirection = -normalize(sunPos-xyz);
  
  //Vector3 sunDirection = normalize(sunPos-xyz);
  Eigen::Vector3f sunDirection = (sunPos-xyz).normalized();
  //float mu_0 = dot_prod(sunDirection,normal);
  float mu_0 = sunDirection.dot(normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  //viewer coordinates relative to the xyz point on the Moon surface
  //Eigen::Vector3f viewDirection = normalize(viewPos-xyz);
  Eigen::Vector3f viewDirection = (viewPos-xyz).normalized();
  //float mu = dot_prod(viewDirection,normal);
  float mu = viewDirection.dot(normal);

  //compute the phase angle /alpha between the viewing direction and the light source direction
  float rad_alpha, deg_alpha;
  float cos_alpha;

  //cos_alpha = dot_prod(sunDirection,viewDirection);
  cos_alpha = sunDirection.dot(viewDirection);
  if (cos_alpha > 1){cos_alpha= 1;}
  if (cos_alpha <-1){cos_alpha=-1;}
  /*
  if ((cos_alpha > 1)||(cos_alpha< -1)){
    cout<<"sunDirection="<<sunDirection<<", viewDirection="<<viewDirection<<endl;
    cout<<"ReflectanceModels.cc::cos_alpha error"<<endl;
  }
  */
  rad_alpha = acos(cos_alpha);
  deg_alpha = rad_alpha*180.0/M_PI;
  //cout<<"ReflectanceModels.cc:: ComputeLunarLambertianReflectanceFromNormal(): deg_alpha="<<deg_alpha<<endl;
  //printf("deg_alpha = %f\n", deg_alpha);

  //Bob Gaskell's model
  //L = exp(-deg_alpha/60.0);

#if 0 // trey
  // perfectly valid for alpha to be greater than 90?
  if (deg_alpha > 90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
  if (deg_alpha < -90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
#endif

  //Alfred McEwen's model
  float A = -0.019;
  float B =  0.000242;//0.242*1e-3;
  float C = -0.00000146;//-1.46*1e-6;

  L = 1.0 + A*deg_alpha + B*deg_alpha*deg_alpha + C*deg_alpha*deg_alpha*deg_alpha;

  //cout<<"ReflectanceModels.cc:: ComputeLunarLambertianReflectanceFromNormal(): L="<<L<<endl;
  
  //        std::cout << " sun direction " << sunDirection << " view direction " << viewDirection << " normal " << normal;
  //        std::cout << " cos_alpha " << cos_alpha << " incident " << mu_0 << " emission " << mu;
  //printf(" deg_alpha = %f, L = %f\n", deg_alpha, L);

  //if (mu_0 < 0.15){ //incidence angle is close to 90 deg
  if (mu_0 < 0.0){
    //mu_0 = 0.15;
    //cout<<"mu_0="<<mu_0<<endl;
    return (0.0);
  }

  if (mu < 0.0){ //emission angle is > 90
      mu = 0.0;
    //return (0.0);
  }

  if (mu_0 + mu == 0){
    //cout<<"ReflectanceModels.cc:: ComputeLunarLambertianReflectanceFromNormal(): Warning division by zero"<<endl;;
    reflectance = 0.0;
  }
  else{
    reflectance = 2*L*mu_0/(mu_0+mu) + (1-L)*mu_0;
    //cout<<"reflectance="<<reflectance<<"mu_0="<<mu_0<<", mu="<<mu<<endl;
  }
  if (reflectance < 0){
    cout<<"ReflectanceModels.cc:: ComputeLunarLambertianReflectanceFromNormal():negative reflectance"<<endl;
    reflectance = 0;
  }
  //cout<<"ReflectanceModels.cc:: ComputeLunarLambertianReflectanceFromNormal(): reflectance="<<reflectance<<endl;
  return reflectance;
}
 
  /*
  float SFS::ComputeLunarLambertianReflectanceFromNormal_1(Eigen::Vector3f const& sunPos, Eigen::Vector3f const& viewPos, 
							   Eigen::Vector3f const& xyz, Eigen::Vector3f const& normal) 
  {

    //float reflectance = ComputeLunarLambertianReflectanceFromNormal(sunPos, viewPos, xyz, normal);
    
    float reflectance;
    float L;
    float alpha;
    
    double len = normal.dot(normal);
    if (abs(len - 1.0) > 1.0e-4){
      std::cerr << "Error: Expecting unit normal in the reflectance computation, in "
		<< __FILE__ << " at line " << __LINE__ << std::endl;
      exit(1);
    }
    
    //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
    //sun coordinates relative to the xyz point on the Moon surface
    //Vector3 sunDirection = normalize(sunPos-xyz);
    Eigen::Vector3f sunDirection = (sunPos-xyz).normalized();
    
  
    float mu_0 = sunDirection.dot(normal);//dot_prod(sunDirection,normal);

    
    //cout<<"mu_0="<<mu_0<<endl;
    double tol = 0.3;
    if (mu_0 < tol){
      // Sun is too low, reflectance is too close to 0, the albedo will be inaccurate
      //cout<<"Low sun angle"<<endl;
      return 0.0;
    }
    
    //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
    //viewer coordinates relative to the xyz point on the Moon surface
    //Vector3 viewDirection = normalize(viewPos-xyz);
    Eigen::Vector3f viewDirection = (viewPos-xyz).normalized();
    float mu = viewDirection.dot(normal);//dot_prod(viewDirection,normal);


    //compute the phase angle (alpha) between the viewing direction and the light source direction
    float deg_alpha;
    float cos_alpha;
    
    cos_alpha = sunDirection.dot(viewDirection);//dot_prod(sunDirection,viewDirection);
    if ((cos_alpha > 1)||(cos_alpha< -1)){
      printf("cos_alpha error\n");
    }
    
    alpha     = acos(cos_alpha);  // phase angle in radians
    deg_alpha = alpha*180.0/M_PI; // phase angle in degrees
    
    //printf("deg_alpha = %f\n", deg_alpha);
    
    //Bob Gaskell's model
    //L = exp(-deg_alpha/60.0);
    
    //Alfred McEwen's model
    float A = -0.019;
    float B =  0.000242;//0.242*1e-3;
    float C = -0.00000146;//-1.46*1e-6;
    
    L = 1.0 + A*deg_alpha + B*deg_alpha*deg_alpha + C*deg_alpha*deg_alpha*deg_alpha;
    
    //printf(" deg_alpha = %f, L = %f\n", deg_alpha, L);
    
    if (mu_0 < 0.0){
      cout<<"mu_0 is negative"<<endl;
      return 0.0;
    }
    
    if (mu < 0.0){ //emission angle is > 90
      //cout<<"Emission angle is > 90"<<endl;
      mu = 0.0;
      //return 0.0;//added by Ara 08/24
    }
    
    if (mu_0 + mu == 0){
      cout<<"Zero sum"<<endl;
      return 0.0;
    }
    else{
      reflectance = 2*L*mu_0/(mu_0+mu) + (1-L)*mu_0;
    }
    if (reflectance <= 0){
      cout<<"Negative reflectance"<<endl;
      return 0.0;
    }
    
    // Attempt to compensate for points on the terrain being too bright
    // if the sun is behind the spacecraft as seen from those points.
    
    //reflectance *= std::max(0.4, exp(-alpha*alpha));
    //reflectance *= ( exp(-phaseCoeffC1*alpha) + phaseCoeffC2 );
    
    return reflectance;
  }
  */
}
