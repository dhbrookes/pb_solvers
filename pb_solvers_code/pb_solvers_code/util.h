//
//  util.hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/24/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef util_hpp
#define util_hpp

#include <stdio.h>
#include "MyMatrix.h"
#include <math.h>
#include <complex>


template <typename T>
class SphPoint;

template<typename T>
class EuPoint;

//Useful typedefs:
typedef EuPoint<double> EPt;
typedef SphPoint<double> ShPt;
typedef complex<double> cmplx;


/*
Right scalar multiplication of a matrix
 */
template<typename T>
const MyMatrix<T> operator*(const MyMatrix<T> mat, const T& rhs)
{
  MyMatrix<T> result = MyMatrix<T>(mat.get_nrows(), mat.get_ncols());
  int i, j;
  for (i = 0; i < mat.get_nrows(); i++)
  {
    for (j= 0; j < mat.get_ncols(); j++)
    {
        result.set_val(i, j, rhs * mat(i, j));
    }
  }
  return result;
}


/*
 lhs scalar multiplication of a matrix
 */
template<typename T>
const MyMatrix<T> operator*(const T& lhs, const MyMatrix<T> mat)
{
  MyMatrix<T> result = MyMatrix<T>(mat.get_nrows(), mat.get_ncols());
  int i, j;
  for (i = 0; i < mat.get_nrows(); i++)
  {
    for (j= 0; j < mat.get_ncols(); j++)
    {
        result.set_val(i, j, lhs * mat(i, j));
    }
  }
  return result;
}



/*
 Class for 3 dimensional euclidean points
 */
template<typename T>
class EuPoint
{
protected:
  T x_;
  T y_;
  T z_;

public:
  EuPoint(T x=T(), T y=T(), T z=T())
  :x_(x), y_(y), z_(z)
  {
  }

  const SphPoint<T> convert_to_spherical() const
  {
    T r, theta, phi;
    r = sqrt(x_*x_ + y_*y_ + z_*z_);
    theta = atan(y_ / x_);
    phi = atan(sqrt(x_*x_ + y_*y_) / z_);
    return SphPoint<T>(r, theta, phi);
  }

  const T& get_x() const  { return x_; }
  const T& get_y() const  { return y_; }
  const T& get_z() const  { return z_; }
};


/*
 Class for 3 dimensional spherical points
 */
template<typename T>
class SphPoint
{
protected:
  T r_;
  T theta_;
  T phi_;
    
public:
  SphPoint(T r=T(), T theta=T(), T phi=T())
  :r_(r), theta_(theta), phi_(phi)
  {
  }
  
  const EuPoint<T> convert_to_euclidean() const
  {
    T x, y, z;
    x = r_ * cos(theta_) * sin(phi_);
    y = r_ * sin(theta_) * sin(phi_);
    z = r_ * cos(phi_);
    return EuPoint<T>(x, y, z);
  }
  
  const T& get_r() const      { return r_;        }
  const T& get_theta() const  { return theta_;    }
  const T& get_phi() const    { return phi_;      }
  
};


#endif /* util_hpp */
