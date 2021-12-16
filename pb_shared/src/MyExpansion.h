//
//  MyExpansion.h
//  pbam
//
//  Created by Lisa Felberg on 8/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef MyExpansion_h
#define MyExpansion_h

#include <map>

namespace pbsolvers
{
  
class ExpansionAccessException: public exception
{
protected:
  int i_, j_;
  int poles_;
  mutable string ss_str_;
  
public:
  ExpansionAccessException(const int i, const int j,
                        const int poles)
  :i_(i), j_(j), poles_(poles), ss_str_()
  {
  }
  
  virtual const char* what() const throw()
  {
    ostringstream ss;
    ss << "Cannot access point [" << i_ << "," <<  j_ <<
    "] in expansion of size " << poles_ << endl;
    ss_str_ = ss.str();
    return ss_str_.c_str();
  }
};

    
class ExpansionArithmeticException: public exception
{
protected:
  int poles1_, poles2_;
  ArithmeticType type_;
  mutable string ss_str_;
      
public:
      
  ExpansionArithmeticException(ArithmeticType type, const int poles1,
                               const int poles2)
  :poles1_(poles1), poles2_(poles2), type_(type), ss_str_()
  {
  }
      
  virtual const char* what() const throw()
  {
    ostringstream ss;
    string start;
    if (type_ == ADDITION) start = "Cannot add matrices of sizes (";
    else if (type_ == MULTIPLICATION)
      start = "Cannot multiply matrices of size (";
    else if (type_ == INNER_PRODUCT)
      start = "Cannot find inner product of vectors of sizes (";
    else start = "Unknown arithmetic error with matrices of sizes (";
    ss << start << poles1_ << " and " << poles2_ << ")" << endl;
    ss_str_ = ss.str();
    return ss_str_.c_str();
  }
};

enum CmplxType { RL, IM};
   
class MyExpansion
{
protected:
  int                 poles_;
  vector<double >      vals_;
  
  map<vector<int>, int> cmplxToDouble_;
  
public:
  
  /*
   Initialize an empty matrix (call default constructors of T)
   Default is 1 x 1 matrix
   */
  MyExpansion(const int poles=1)
  : poles_(poles), vals_ (poles*poles, 0.0)
  {
    map_cmplx_to_double();
  }
  
  MyExpansion(vector<double> vals)
  : vals_(vals), poles_( (int) sqrt((int) vals.size()))
  {
    map_cmplx_to_double();
  }
  
  MyExpansion(const int poles, double * val)
  :poles_(poles), vals_(poles)
  {
    int i;
    for (i = 0; i < poles*poles; i++)
      vals_[i] = val[i];
    map_cmplx_to_double();
  }
  /*
   Fill with a default value (good for initializing memory)
   */
  MyExpansion(const int poles, double default_val)
  :poles_(poles), vals_(poles*poles, default_val)
  {
    map_cmplx_to_double();
  }
  
  void map_cmplx_to_double()
  {
    int i, j, ct(0);
    vector<int> cmplxMp(3);
    
    for ( i = 0; i < poles_; i++)
      for( j = 0; j <= i; j++)
      {
        cmplxMp[0] = i; cmplxMp[1] = j;
        cmplxMp[2] = (int) RL;
        cmplxToDouble_[cmplxMp] = ct;
        ct++;
        
        if (j > 0)
        {
          cmplxMp[0] = i; cmplxMp[1] = j;
          cmplxMp[2] = (int)IM;
          cmplxToDouble_[cmplxMp] = ct;
          ct++;
        }
      }
  }
  
  /*
   Set the value of a coordinate in this matrix given the position (i, j).real
   */
  void set_val(const int i, const int j, CmplxType type, double val)
  {
    vector<int> cmplxMp = {i,j,type};
    vals_[cmplxToDouble_[cmplxMp]] = val;
  }
  
  /*
   Set the value of a coordinate in this matrix given the position (i, j).real
   */
  void set_val_cmplx(const int i, const int j, complex<double> val)
  {
    vector<int> cmplxMp = {i,j,RL};
    vals_[cmplxToDouble_[cmplxMp]] = val.real();
    
    if (j > 0) vals_[cmplxToDouble_[cmplxMp]+1] = val.imag();
  }
  
  void set_val(const int i, double val) { vals_[i] = val; }
  
  /*
   Element access operator given position (i, j) and real/imag
   */
  double operator()(const int i, const int j, CmplxType type)
  {
    if (i > poles_ || j > i)
    {
      throw ExpansionAccessException(i, j, poles_);
    }
    else
    {
      
      if ((j==0) && (type==IM))
        return 0.0;
      
      if ((j < 0) && (type==IM))
      {
        vector<int> cmplxMp = {i,-j,type};
        return -vals_[cmplxToDouble_[cmplxMp]];
      } else
      {
        vector<int> cmplxMp = {i,abs(j),type};
        return vals_[cmplxToDouble_[cmplxMp]];
      }
    }
  }
  
  /*
   Element access operator given position (i, j) and real/imag
   */
  double operator()(const int i)
  {
    if (i < 0 || i > poles_*poles_)
    {
      throw ExpansionAccessException(i, i, poles_);
    }
    else
    {
      return vals_[i];
    }
  }
  
  /*
   Element access operator given position (i, j) and real/imag
   */
  complex<double> get_cmplx(const int i, const int j)
  {
    if (i > poles_ || j > i)
    {
      throw ExpansionAccessException(i, j, poles_);
    }
    else
    {
      vector<int> cmplxMp = {i,abs(j),RL};
      int locMap = cmplxToDouble_[cmplxMp];
      if (j==0)
        return complex<double>(vals_[locMap],0.0);
      else if (j<0)
        return complex<double>(vals_[locMap], -vals_[locMap+1]);
      
      return complex<double>(vals_[locMap], vals_[locMap+1]);
    }
  }
  
  
  /*
   Addition operator returns new matrix
   */
  MyExpansion operator+(MyExpansion& rhs)
  {
    if (poles_ != rhs.poles_)
    {
      throw ExpansionArithmeticException(ADDITION, poles_, rhs.poles_);
    }
    
    MyExpansion result = MyExpansion(poles_);
    int i;
    for (i = 0; i < poles_*poles_; i++)
    {
      result.set_val(i, vals_[i] + rhs(i));
    }
    return result;
  }
  
  /*
   summation operator adds to existing matrix
   */
  MyExpansion operator+=(MyExpansion& rhs)
  {
    if (poles_ != rhs.poles_)
    {
      throw ExpansionArithmeticException(ADDITION, poles_, rhs.poles_);
    }
    int i;
    for (i = 0; i < poles_*poles_; i++)
    {
      set_val(i, vals_[i] + rhs(i));
    }
    return *this;
  }
  
  void print_expansion(int p)
  {
    int i, j;
    for (i = 0; i < p; i++)
    {
      for (j = 0; j <= i ; j++)
      {
        double  r = this->operator()(i, j, RL);
        double im = this->operator()(i, j, IM);
        r  = fabs( r) > 1e-15 ?  r : 0;
        im = fabs(im) > 1e-15 ? im : 0;
        cout << "("<< setprecision(9) << r << ", " << im <<") ";
      }
      cout << endl;
    }
  }

  /*
   Vector Expansion multiplication. If this has n poles, then rhs 
   must be a vector of length n.
   */
  MyExpansion operator*(MyVector<double> & rhs)
  {
    if (poles_ != rhs.get_nrows())
      throw ExpansionArithmeticException(MULTIPLICATION, poles_,
                                         (int)rhs.get_nrows());

    MyExpansion result = MyExpansion(poles_);
    int i, j, ct(0);
    for (i = 0; i < poles_; i++)
    {
      for (j= 0; j < 2*i+1; j++) //Accounting for real+imag
      {
        result.set_val(ct, vals_[ct] * rhs[i]);
        ct++;
      }
    }
    return result;
  }
  
  /*
   Vector Expansion multiplication. If this has n poles, then rhs
   must be a vector of length n.
   */
  MyExpansion operator*(double scalar)
  {
    MyExpansion result = MyExpansion(poles_);
    int i, j, ct(0);
    for (i = 0; i < poles_; i++)
    {
      for (j= 0; j < 2*i+1; j++) //Accounting for real+imag
      {
        result.set_val(ct, vals_[ct] * scalar);
        ct++;
      }
    }
    return result;
  }
  
  const int get_poles() { return poles_; }
  
};
    
    
class MyGradExpansion
{
protected:
  int                 poles_;
  vector<MyExpansion>  exps_;
      
  map<vector<int>, int> cmplxToDouble_;
      
public:
  
  MyGradExpansion(int poles = 1) : exps_(3, MyExpansion(poles)), poles_(poles)
  {
  }
  
  MyExpansion get_dim(int dim) { return exps_[dim];}
  
  complex<double> get_dimij_cmplx(int dim, int i, int j)
  { return exps_[dim].get_cmplx(i, j);}
  
  double get_dimi(int dim, int i)
  { return exps_[dim](i);}
  
  void set_dim(int dim, MyExpansion exp) { exps_[dim] = exp;}
  
  void set_dim_val_cmplx(int dim, int i, int j, complex<double> val)
  { exps_[dim].set_val_cmplx(i, j, val);}
  
  void set_dimi(int dim, int i, double val) { exps_[dim].set_val(i, val);}
  
  /*
   Addition operator returns new matrix
   */
  MyGradExpansion operator+(MyGradExpansion &rhs)
  {
    if (poles_ != rhs.poles_)
    {
      throw ExpansionArithmeticException(ADDITION, poles_, rhs.poles_);
    }
    
    MyGradExpansion result = MyGradExpansion(poles_);
    int i,j;
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j<poles_*poles_; j++)
        result.set_dimi(i, j, exps_[i](j) + rhs.exps_[i](j));
    }
    return result;
  }
  
  /*
   summation operator adds to existing matrix
   */
  MyGradExpansion operator+=(MyGradExpansion & rhs)
  {
    if (poles_ != rhs.poles_)
    {
      throw ExpansionArithmeticException(ADDITION, poles_, rhs.poles_);
    }
    int i, j;
    for (i = 0; i < 3; i++)
    {
      for (j = 0; j<poles_*poles_; j++)
      exps_[i].set_val(j, exps_[i](j) + rhs.exps_[i](j));
    }
    return *this;
  }
  
  
};

} /* namespace pbsolvers */
#endif /* MyExpansion_h */
