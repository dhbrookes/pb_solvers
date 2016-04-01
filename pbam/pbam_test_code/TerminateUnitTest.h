//
//  TerminateUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/31/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef TerminateUnitTest_h
#define TerminateUnitTest_h

#include "BD.h"

class TermUTest : public ::testing::Test
{
public :
  
protected :

  
  virtual void SetUp() { }
  
  virtual void TearDown() {}
} ; // end BDUTest


// Checking time termination class
TEST_F(TermUTest, timeTerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }

  auto sys = make_shared<System> ( mol );
  shared_ptr<TimeTerminate> term = make_shared<TimeTerminate>(30);
  
  EXPECT_NEAR(sys->get_time(), 0, preclim);
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->set_time( 40);
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking x<= termination class
TEST_F(TermUTest, xLETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(1, X, LEQ, -10);

  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-5, 0, 0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-5, 0, 0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking x>= termination class
TEST_F(TermUTest, xGETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(2, X, GEQ, 12.54);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(5, 0, 0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(7.3, 0, 0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking y<= termination class
TEST_F(TermUTest, yLETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(0, Y, LEQ, -10);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(-5, 0, 0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(-5, -10, 0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking y>= termination class
TEST_F(TermUTest, yGETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(1, Y, GEQ, 10);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(.5, -2, 0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-5, 7, 0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking z<= termination class
TEST_F(TermUTest, zLETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(0, Z, LEQ, -50);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(100, 0, -5));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(-45, -10, -45));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking z>= termination class
TEST_F(TermUTest, zGETerm)
{
  const int nmol = 3;
  vector<Molecule> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    Molecule molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(2, Z, GEQ, 29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 10.4, 20));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), true);
}

#endif /* TerminateUnitTest_h */
