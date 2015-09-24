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





#endif /* util_hpp */
