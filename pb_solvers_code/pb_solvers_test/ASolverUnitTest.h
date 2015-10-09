//
//  ASolverUnitTest.h
//  pbsolvers
//
//  Created by Marielle Soniat on 10/6/15.
//  Copyright (c) 2015 Marielle Soniat. All rights reserved.
//

#ifndef pbsolvers_ASolverUnitTest_h
#define pbsolvers_ASolverUnitTest_h

#include <iostream>
#include "ASolver.h"
using namespace std;


class ASolverUTest : public ::testing::Test
{
public :
  
protected :
  
  int vals_;
  
  Constants const_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    
    mol_.clear( );
    EPt pos[2] = { EPt( 0.0, 0.0, -5.0 ), EPt( 0.0, 0.0, 5.0 ) };
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      double a = 1.0;
      vector<double> charges(1);
      vector<EPt> posCharges(1);
      
      charges[0] = 1.0;
      posCharges[0] = pos[molInd];
      
      Molecule molNew = Molecule( M, a, charges, posCharges );
      mol_.push_back( molNew );
    }
    
  }
  virtual void TearDown() {}
} ;

// functions to test:
/*
 compute_gamma();
 compute_delta();
 compute_E();
 calc_mol_sh
 calc_indi_gamma
 calc_indi_delta
 calc_indi_e
 */

TEST_F(ASolverUTest, checkVals)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( vals );
  BesselCalc bCalcu        = BesselCalc( vals, &bConsta);
  SHCalcConstants SHConsta = SHCalcConstants( vals );
  SHCalc SHCalcu           = SHCalc( vals, &SHConsta );
  System sys               = System( const_, mol_ );
  ASolver ASolvTest = ASolver( nmol, vals, &bCalcu, &SHCalcu, sys);
  
  for (int molIt = 0; molIt < nmol; molIt++ )
  {
    //ASolvTest.calc_mol_sh( sys.get_molecule( molIt ));
  }
}


/*
 
 // testing matrix multiplication
 TEST_F(MyMatrixUTest, matMul)
 {
 MyMatrix <double> testMat1( testInp_ );
 MyMatrix <double> testMat2( testInp2_ );
 MyMatrix <double> testMat3;
 
 testMat3 = testMat1 * testMat2;
 ASSERT_EQ( testMat3.get_nrows() , 2 );
 ASSERT_EQ( testMat3.get_ncols() , 2 );
 
 EXPECT_NEAR( testMat3(0,0),  50.0, preclim);
 EXPECT_NEAR( testMat3(0,1),  60.0, preclim);
 EXPECT_NEAR( testMat3(1,0), 114.0, preclim);
 EXPECT_NEAR( testMat3(1,1), 140.0, preclim);
 }
 
 // testing matrix assertions for accessing matrix members
 TEST_F(MyMatrixUTest, exceptAccess)
 {
 MyMatrix <double> testMat_( testInp_ );
 ASSERT_THROW(testMat_( 0, 5), MatrixAccessException);
 ASSERT_THROW(testMat_( 5, 0), MatrixAccessException);
 ASSERT_THROW(testMat_(-1, 0), MatrixAccessException);
 ASSERT_THROW(testMat_( 0,-1), MatrixAccessException);
 }
 
 
 */



#endif
