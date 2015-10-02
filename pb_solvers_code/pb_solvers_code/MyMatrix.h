//
//  MyMatrix.h
//  pb_solvers_code
//
//  Created by David Brookes on 9/24/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef MyMatrix_h
#define MyMatrix_h

#include <stdio.h>
#include <vector>
#include <sstream>


using namespace std;


class MatrixAccessException: public exception
{
protected:
    int i_, j_;
    int nrows_, ncols_;
    
public:
    MatrixAccessException(const int i, const int j, const int nrows, const int ncols)
    :i_(i), j_(j), nrows_(nrows), ncols_(ncols)
    {
    }
    
    virtual const char* what() const throw()
    {
        ostringstream ss;
        ss << "Cannot access point [" << i_ << "," <<  j_ << "] in matrix of size (" << nrows_ << "," << ncols_ << ")" << endl;
        return ss.str().c_str();
    }
};


enum ArithmeticType { ADDITION, MULTIPLICATION, INNER_PRODUCT };
        
        
class MatrixArithmeticException: public exception
{
protected:
    int nrows1_, ncols1_;
    int nrows2_, ncols2_;
    ArithmeticType type_;
    
public:
    
    MatrixArithmeticException(ArithmeticType type, const int nrows1, const int ncols1,
                              const int nrows2, const int ncols2)
    :nrows1_(nrows1), ncols1_(ncols1), nrows2_(nrows2), ncols2_(ncols2), type_(type)
    {
    }
    
    virtual const char* what() const throw()
    {
        ostringstream ss;
        string start;
        if (type_ == ADDITION) start = "Cannot add matrices of sizes (";
        else if (type_ == MULTIPLICATION) start = "Cannot multiply matrices of size (";
        else if (type_ == INNER_PRODUCT) start = "Cannot find inner product of vectors of sizes (";
        else start = "Unknown arithmetic error with matrices of sizes (";
        ss << start << nrows1_ << "," << ncols1_ << ") and (" << nrows2_ << "," << ncols2_ << ")" << endl;
        return ss.str().c_str();
    }

};

        
template <typename T>
class MyMatrix
{
protected:
    int                 nrows_;
    int                 ncols_;
    vector< vector<T> >  vals_;  //length of first vector is number of rows, length of second is ncols
    

public:

    /*
     Initialize an empty matrix (call default constructors of T)
     Default is 1 x 1 matrix
     */
    MyMatrix(const int nrows=1, const int ncols=1)
    : nrows_(nrows), ncols_(ncols), vals_ (nrows, vector<T>(ncols))
    {
    }
    
    MyMatrix(vector< vector<T> >& vals)
    : vals_(vals), nrows_(vals.size()), ncols_(vals[0].size())
    {
    }
    
    /*
     Set the value of a coordinate in this matrix given the position (i, j)
     */
    void set_val(const int i, const int j, const T& val)
    {
        vals_[i][j] = val;
    }
    
    /*
     Element access operator given position (i, j)
     */
    const T& operator()(const int i, const int j) const
    {
        if (i < 0 || j < 0 || i > nrows_ || j > ncols_)
        {
            throw MatrixAccessException(i, j, nrows_, ncols_);
        }
        else
        {
            return vals_[i][j];
        }
    }
    
    /*
     Addition operator returns new matrix
     */
    const MyMatrix<T> operator+(const MyMatrix<T>& rhs) const
    {
        if (ncols_ != rhs.ncols_ || nrows_ != rhs.nrows_)
        {
            throw MatrixArithmeticException(ADDITION, nrows_, ncols_, rhs.nrows_, rhs.ncols_);
        }
        
        MyMatrix<T> result = MyMatrix<T>(nrows_, ncols_);
        int i, j;
        for (i = 0; i < nrows_; i++)
        {
            for (j= 0; j < ncols_; j++)
            {
                result.set_val(i, j, this(i, j) + rhs(i, j));
            }
        }
        return result;
    }
    
    /*
     Matrix multiplication. If this is size n x m, then rhs must be size m x p
     */
    const MyMatrix<T> operator*(const MyMatrix<T>& rhs) const
    {
        if (ncols_ != rhs.nrows_)
        {
            throw MatrixArithmeticException(MULTIPLICATION, nrows_, ncols_, rhs.nrows_, rhs.ncols_);
        }
        
        int n, m, p;
        n = nrows_;
        m = ncols_;
        p = rhs.ncols_;
        
        MyMatrix<T> result = MyMatrix<T>(n, p);
        int i, j, k;
        T inner_sum;
        for (i = 0; i < n; i++)
        {
            for (j= 0; j < p; j++)
            {
                inner_sum = T();  // default constructor should be equivalent to zero
                for (k = 0; k < m; k++)
                {
                    inner_sum += this(i, k) * rhs(k, j);
                }
                result.set_val(i, j, inner_sum);
            }
        }
        return result;
    }
    
        
    const int get_nrows() const { return nrows_; }
    const int get_ncols() const { return ncols_; }

};
        
        
/*
 Vector class is implemented as extension of matrix class but with only one column
 */
template<typename T>
class MyVector : public MyMatrix<T>
{
public:
    /*
     Initialize empty vector given the size
     */
    MyVector(const int size=1)
    :MyMatrix<T>(size, 1)
    {
    }
    
    MyVector(const vector<T>& vals)
    :MyMatrix<T>(vals.size(), 1)
    {
        int i;
        for (i = 0; i < vals.size(); i++)
        {
            this->vals_[i][0] = vals[0];
        }
    }
    
    void set_val(const int i, const T& val)
    {
        this->MyMatrix<T>::set_val(i, 0, val);
    }
    
    /*
     Access operator with brackets only requires one value
     */
    const T& operator[](int i)
    {
        return this(i, 0);
    }
    
    /*
     The multiplication operator now computes the inner product
     */
    const T operator*(const MyVector<T>& rhs)
    {
        if (rhs->nrows_ != this->nrows)
        {
            throw MatrixArithmeticException(INNER_PRODUCT, this->nrows_, this->ncols_,
                                            rhs->nrows_, rhs->ncols_);
        }
        T out = T();
        int i;
        for(i = 0; i < this->nrows_; i++)
        {
            out = out + (this[i] * rhs[i]);
        }
        return out;
    }
};

        
/*
 Convenience classes for matrices of matrices, vectors of vectors
 and vectors of matrices.
 
 These are used by calling the name of the class and then ::type
 to retrieve the typedef that they contain.
 
 So to get a matrix of matrices containing ints you would write:
 
 MatofMats<int>::type my_mat;
 
 */
template <typename T>
struct MatOfMats
{
    typedef MyMatrix<MyMatrix<T> > type;
};

template <typename T>
struct VecOfVecs
{
    typedef MyVector<MyVector<T> > type;
};
        
template <typename T>
struct VecOfMats
{
    typedef MyVector<MyMatrix<T> > type;
};


        
#endif /* MyMatrix_h */
        
