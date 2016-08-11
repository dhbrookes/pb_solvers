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


class ExpansionAccessException: public exception
{
protected:
  int i_, j_;
  int poles_;
  
public:
  ExpansionAccessException(const int i, const int j,
                        const int poles)
  :i_(i), j_(j), poles_(poles)
  {
  }
  
  virtual const char* what() const throw()
  {
    ostringstream ss;
    ss << "Cannot access point [" << i_ << "," <<  j_ <<
    "] in expansion of size " << poles_ << endl;
    return ss.str().c_str();
  }
};

    
class ExpansionArithmeticException: public exception
{
protected:
  int poles1_, poles2_;
  ArithmeticType type_;
      
public:
      
  ExpansionArithmeticException(ArithmeticType type, const int poles1,
                               const int poles2)
  :poles1_(poles1), poles2_(poles2), type_(type)
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
    return ss.str().c_str();
  }
};

enum CmplxType { REAL, IMAG };
   
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
  : poles_(poles), vals_ (poles*poles)
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
        cmplxMp = {i,j,REAL};
        cmplxToDouble_[cmplxMp] = ct;
        ct++;
        
        if (j > 0)
        {
          cmplxMp = {i,j,IMAG};
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
    vector<int> cmplxMp = {i,j,REAL};
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
      
      if ((j==0) && (type==IMAG))
        return 0.0;
      
      if ((j < 0) && (type==IMAG))
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
      vector<int> cmplxMp = {i,abs(j),REAL};
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
  MyExpansion operator+(MyExpansion rhs)
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
  MyExpansion operator+=(MyExpansion rhs)
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
        double  r = this->operator()(i, j, REAL);
        double im = this->operator()(i, j, IMAG);
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
  
  const int get_poles() { return poles_; }
  
};
    
    
class MyGradExpansion
{
protected:
  int                 poles_;
  vector<MyExpansion>  exps_;
      
  map<vector<int>, int> cmplxToDouble_;
      
public:
  
  MyGradExpansion(int poles) : exps_(3, MyExpansion(poles)), poles_(poles)
  {
  }
  
  MyExpansion get_dim(int dim) { return exps_[dim];} 
  
  
  
  
};

#endif /* MyExpansion_h */
