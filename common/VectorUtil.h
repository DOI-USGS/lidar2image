#ifndef ATK_COMMON_VECTOR_UTIL_H
#define ATK_COMMON_VECTOR_UTIL_H

#include <Eigen/Dense>

namespace at
{

  //Eigen 3.2+ defines .hasNaN() and .allFinite() (briefly called isFinite), but
  //until we're using at least that version everywhere, we need our own
  //To keep life simple, these are defined generically but will only compile
  //of course on Eigen types
  template<class VecType>
    inline bool VectorHasNaN(const VecType &vec)
    {
      //cute parallel trick copied from Eigen source: NaN doesn't equal anything
      return !((vec.array()==vec.array()).all());
    }

  template<class VecType>
    inline bool VectorIsFinite(const VecType &vec)
    {
      //cute parallel trick copied from Eigen source: inf-inf = nan, nan-nan=nan
      return !(VectorHasNaN(vec-vec));
    }

};

#endif /* ATK_COMMON_VECTOR_UTIL_H */

