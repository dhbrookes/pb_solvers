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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
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

// Checking z<= termination class with many MoleculeAM of one type
TEST_F(TermUTest, zLETermManyType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  vector<int> typ = { 0, 0, 1}; vector<int> typind = { 0, 1, 0};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(0, Z, LEQ, -50);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(100, 0, -500));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-45, -10, -55));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking z>= termination class for many instances of one type
TEST_F(TermUTest, xGETermManyType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(1, X, GEQ, 29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 10.4, 20));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(2, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking x<= class for many instances of one type with PBCs
TEST_F(TermUTest, xLETermManyTypePBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  auto term = make_shared<CoordTerminate>(0, X, LEQ, -29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-4.5, 9.5, 0));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-10, -10.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(-38, -44, 0.0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking y>= class for many instances of one type with PBCs
TEST_F(TermUTest, yGETermManyTypePBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  auto term = make_shared<CoordTerminate>(1, Y, GEQ, 39.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-4.5, 9.5, 0));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(-21, -10.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(0, Pt(28, 24, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(18, 44, 0.0));
  EXPECT_EQ( term->is_terminated( sys), true);
}


// Checking r termination class
TEST_F(TermUTest, rTermType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(1, R, GEQ, 29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 10.4, 20));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking r termination class for many instances of one type
TEST_F(TermUTest, rTermManyType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  auto term = make_shared<CoordTerminate>(1, R, GEQ, 29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 1.4, 5));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(2, Pt(35, 0.5, 9.7));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking r class for many instances of one type with PBCs
TEST_F(TermUTest, rTermManyTypePBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  auto term = make_shared<CoordTerminate>(1, R, GEQ, 29.4);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-4.5, 9.5, 0));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(38, -44, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(2, Pt(-4.5, 9.5, 0));
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-14.5, 9.5, 0));
  EXPECT_EQ( term->is_terminated( sys), true);
}


// Checking contact termination class
TEST_F(TermUTest, contTermType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    molInd, 0);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  vector<int> cont_pair = { 1, 2};
  auto term = make_shared<ContactTerminate>(cont_pair, 1.0);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(0, Pt(-5, 10.4, 20));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-5, 0.0, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(0, -2.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking contact class for many instances of one type
TEST_F(TermUTest, contTermManyType)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol );
  vector<int> cont_pair = { 0, 1};
  auto term = make_shared<ContactTerminate>(cont_pair, 1.0);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-2.5, 0, 0));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(0, -0.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(0, -1.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), true);
}

// Checking contact class for many instances of one type with PBCs
TEST_F(TermUTest, contTermManyTypePBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  vector<int> cont_pair = { 0, 1};
  auto term = make_shared<ContactTerminate>(cont_pair, 1.0);
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(2, Pt(-24.5, 19.5, 0));
  
  EXPECT_EQ( term->is_terminated( sys), false);
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  EXPECT_EQ( term->is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( term->is_terminated( sys), true);

}


// Checking contact + time class combined with OR
TEST_F(TermUTest, combineTermORPBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  vector<int> cont_pair = { 0, 1};
  vector<shared_ptr<BaseTerminate > > all_term;
  all_term.push_back(make_shared<ContactTerminate>(cont_pair, 1.0));
  all_term.push_back(make_shared<TimeTerminate>(100.0));
  
  CombineTerminate combine(all_term, ONE);
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(2, Pt(-24.5, 19.5, 0));
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), true);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 40);
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 102.1);
  EXPECT_EQ( combine.is_terminated( sys), true);
}


// Checking contact + time class combined with AND
TEST_F(TermUTest, combineTermANDPBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  vector<int> cont_pair = { 0, 1};
  vector<shared_ptr<BaseTerminate > > all_term;
  all_term.push_back(make_shared<ContactTerminate>(cont_pair, 1.0));
  all_term.push_back(make_shared<TimeTerminate>(100.0));
  
  CombineTerminate combine(all_term, ALL);
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(2, Pt(-24.5, 19.5, 0));
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 40);
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 102.1);
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(-8, 14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), true);
}

// Checking x<= + contact + time class combined with AND
TEST_F(TermUTest, combineTerm3ORPBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  vector<int> cont_pair = { 0, 1};
  vector<shared_ptr<BaseTerminate > > all_term;
  all_term.push_back(make_shared<ContactTerminate>(cont_pair, 1.0));
  all_term.push_back(make_shared<CoordTerminate>(1, X, LEQ, -20.4));
  all_term.push_back(make_shared<TimeTerminate>(100.0));
  
  CombineTerminate combine(all_term, ONE);
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(0, Pt(-24.5, 19.5, 0));
  EXPECT_EQ( combine.is_terminated( sys), true);
  
  sys->translate_mol(0, Pt( 24.5, 19.5, 0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  sys->translate_mol(2, Pt(-25, 19.5, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), true);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 102.1);
  EXPECT_EQ( combine.is_terminated( sys), true);
  
  sys->set_time( 40);
  EXPECT_EQ( combine.is_terminated( sys), false);
}

// Checking x<= + contact + time class combined with AND
TEST_F(TermUTest, combineTerm3ANDPBC)
{
  const int nmol = 3;
  vector<MoleculeAM> mol;
  vector<int> typ = { 1, 0, 1}; vector<int> typind = { 0, 0, 1};
  Pt pos[nmol] = {Pt(0.0, 0.0, 0.0), Pt(0.0, 5.0, 0.0), Pt(5, 0, 0)};
  for (int molInd = 0; molInd < nmol; molInd ++ )
  {
    int M = 1;
    vector<double> charges(M); vector<double> vdW(M); vector<Pt> posCharges(M);
    charges[0] = 2.0; posCharges[0] = pos[molInd]; vdW[0] = 0.0;
    
    MoleculeAM molNew( "stat", 1.0, charges, posCharges, vdW, pos[molInd],
                    typ[molInd], typind[molInd]);
    mol.push_back( molNew );
  }
  
  auto sys = make_shared<System> ( mol, 20, 40 );
  vector<int> cont_pair = { 0, 1};
  vector<shared_ptr<BaseTerminate > > all_term;
  all_term.push_back(make_shared<ContactTerminate>(cont_pair, 1.0));
  all_term.push_back(make_shared<CoordTerminate>(1, X, LEQ, -20.4));
  all_term.push_back(make_shared<TimeTerminate>(100.0));
  
  CombineTerminate combine(all_term, ALL);
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(0, Pt(-24.5, 19.5, 0));
  EXPECT_EQ( all_term[1]->is_terminated( sys), true);
  
  EXPECT_EQ( combine.is_terminated( sys), false);
  sys->translate_mol(1, Pt(10, -10.5, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(8, -14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 40);
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->translate_mol(1, Pt(-8, 14, 0.0));
  EXPECT_EQ( combine.is_terminated( sys), false);
  
  sys->set_time( 102.1);
  EXPECT_EQ( combine.is_terminated( sys), true);
}


#endif /* TerminateUnitTest_h */
