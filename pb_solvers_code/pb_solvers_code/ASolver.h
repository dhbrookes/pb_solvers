//
//  ASolver.hpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/25/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef ASolver_h
#define ASolver_h

#include <stdio.h>
#include "MyMatrix.h"


/*
 This class is designed to compute the vector A defined on pg. 544 of Lotan, Head-Gordon 2006
 */
class ASolver
{
protected:
    typedef MyMatrix<MyMatrix<double> > MatOfMats;
    typedef MyVector<MyVector<double> > VecOfVecs;
    
    VecOfVecs A_, E_;
    MatOfMats gamma_, delta_, T;
    
    
};

#endif /* ASolver_hpp */
