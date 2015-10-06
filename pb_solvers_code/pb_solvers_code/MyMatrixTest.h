//
//  MyMatrixTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/6/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#ifndef MyMatrixTest_h
#define MyMatrixTest_h

#include "MyMatrix.h"
using namespace std;

/*
 Class for testing matrix class
 */
class MyMatrixTest
{
public:
  void TestMatrix( )
  {
    const int nCol = 4;
    const int nRow = 2;
    //double precLim = 1.0e-4;
    
    double val = 1.0;
    vector< vector<double> > testInp;
    testInp.resize(nRow);
    
    for ( int i = 0.0; i < nRow; i++ )
    {
      testInp[i].resize(nCol);
      for ( int j = 0.0; j < nCol; j++ )
      {
        testInp[i][j] = val;
        val += 1.0;
      }
    }
    
    MyMatrix <double> Mat( testInp );
    
//    BesselConstants bConstTest = BesselConstants( nPol );
//    assert( bConstTest.get_n() == nPol ); // make sure numVals stores right
//
//    
//    for (int i = 0; i < nPol; i++) // check that our prefactors are right
//      assert( abs(kPreFactors[i]-bConstTest.get_kconst_val( i )) < precLim);
  }
}; // end MyMatrixTest

#endif /* MyMatrixTest_h */
