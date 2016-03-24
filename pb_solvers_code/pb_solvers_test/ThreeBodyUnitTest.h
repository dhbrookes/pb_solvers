//
//  ThreeBodyUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef ThreeBodyUnitTest_h
#define ThreeBodyUnitTest_h

#include "ThreeBody.h"

class TBDUTest : public ::testing::Test
{
  public :
  
  protected :
  
  int vals_;
  vector< Molecule > mol_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -8.0 ), Pt( 0, 0, 0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -7.0 ), Pt( 0, 0, 1 ) };
    double cg[2] = { -5.0, -5.0};
    double rd[2] = {  3.6,  3.6};
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 3; vector<double> charges(M);
      vector<double> vdW(M); vector<Pt> posCharges(M);
      charges[0]=cg[molInd]; vdW[0]=0.0; posCharges[0]=cgPos[molInd];
      charges[1]=cg[molInd]; vdW[1]=0.0; posCharges[1]=pos[molInd]+Pt(1,0,0);
      charges[2]=cg[molInd]; vdW[2]=0.0; posCharges[2]=pos[molInd]+Pt(0,1,0);
      
      Molecule molNew("stat",rd[molInd],charges,posCharges,vdW,pos[molInd]);
      mol_.push_back( molNew );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  
};

TEST_F(TBDUTest, computeGroups)
{
  vector<Molecule> mol_sing_;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < 9; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol_sing_.push_back( molNew );
  }
  
  const int vals = 18;
  Constants const_( INTERNAL );
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol_sing_);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  EXPECT_TRUE( dim.size() == 36 );
  EXPECT_TRUE( tri.size() == 84 );
  
  EXPECT_TRUE( dim[0][0] == 0 );
  EXPECT_TRUE( dim[0][1] == 1 );
  
  EXPECT_TRUE( dim[35][0] == 7 );
  EXPECT_TRUE( dim[35][1] == 8 );
  
  EXPECT_TRUE( tri[0][0] == 0 );
  EXPECT_TRUE( tri[0][1] == 1 );
  EXPECT_TRUE( tri[0][2] == 2 );
  
  EXPECT_TRUE( tri[83][0] == 6 );
  EXPECT_TRUE( tri[83][1] == 7 );
  EXPECT_TRUE( tri[83][2] == 8 );

}

TEST_F(TBDUTest, twoBD)
{
  vector<Molecule> mol;
  const int num = 7;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  threeBodTest.solveNmer(2, "");
  threeBodTest.calcTwoBDEnForTor();
  
  int j;
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << threeBodTest.get_energyi_approx(j) << "\t"
    << threeBodTest.get_forcei_approx(j).x();
    cout << ", " << threeBodTest.get_forcei_approx(j).y()<< ", ";
    cout << threeBodTest.get_forcei_approx(j).z()<< "\t"<< endl;
  }
  
  cout << endl;
  
  shared_ptr<ASolver> aSolvall = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  aSolvall->solve_A(1e-4); aSolvall->solve_gradA(1e-4);
  PhysCalc calcphys( aSolvall);
  calcphys.calc_all();
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << calcphys.get_omegai(j) << "\t"
    << calcphys.get_forcei(j)[0];
    cout << ", " << calcphys.get_forcei(j)[1]<< ", ";
    cout << calcphys.get_forcei(j)[2]<< "\t" << endl;
  }
  
  cout << endl;
  
}


TEST_F(TBDUTest, threeBD)
{
  vector<Molecule> mol;
  const int num = 4;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  threeBodTest.solveNmer(2, "");
  threeBodTest.solveNmer(3, "");
  threeBodTest.calcTBDEnForTor();
  
  int j;
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << threeBodTest.get_energyi_approx(j) << "\t"
    << threeBodTest.get_forcei_approx(j).x();
    cout << ", " << threeBodTest.get_forcei_approx(j).y()<< ", ";
    cout << threeBodTest.get_forcei_approx(j).z()<< "\t"<< endl;
  }
  
  cout << endl;
  
  shared_ptr<ASolver> aSolvall = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                       make_shared<Constants>
                                                       (const_), vals);
  aSolvall->solve_A(1e-4); aSolvall->solve_gradA(1e-4);
  PhysCalc calcphys( aSolvall);
  calcphys.calc_all();
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << calcphys.get_omegai(j) << "\t"
    << calcphys.get_forcei(j)[0];
    cout << ", " << calcphys.get_forcei(j)[1]<< ", ";
    cout << calcphys.get_forcei(j)[2]<< "\t" << endl;
  }
  
  cout << endl;
}


TEST_F(TBDUTest, threeBDfor3)
{
  vector<Molecule> mol;
  const int num = 3;
  Pt pos[9] = {  Pt( 0.0, 0.0, 0.0 ),Pt( 5.0, 0.0, 0.0 ),Pt( -5.0, 0.0, 0.0 ),
    Pt( -5.0, -5.0, 0.0 ),Pt( -5.0, 5.0, 0.0),Pt( 5.0, -5.0, 0.0 ),
    Pt( 5.0, 5.0, 0.0 ),Pt( 0.0, -5.0, 0.0),Pt( 0.0, 5.0, 0.0),};
  for (int molInd = 0; molInd < num; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<double> vdW(M);
    vector<Pt> posCharges(M);
    charges[0]=2.0;  vdW[0]=0;posCharges[0]=pos[molInd];
    charges[1]=-2.0; vdW[1]=0;posCharges[1]=pos[molInd]+Pt(1.0, 0.0, 0.0);
    charges[2]=2.0;  vdW[2]=0;posCharges[2]=pos[molInd]+Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", 2.0, charges, posCharges, vdW, pos[molInd]);
    mol.push_back( molNew );
  }
  
  const int vals = 10;
  Constants const_( INTERNAL );
  const_.set_salt_concentration(0.1);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*vals);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*vals, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*vals);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*vals, SHConsta);
  shared_ptr<System> sys = make_shared<System>(mol);
  shared_ptr<ASolver> ASolvTest = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        make_shared<Constants>
                                                        (const_), vals);
  ThreeBody threeBodTest( ASolvTest );
  vector<vector<int > > dim = threeBodTest.getDimers();
  vector<vector<int > > tri = threeBodTest.getTrimers();
  
  threeBodTest.solveNmer(2, "");
  threeBodTest.solveNmer(3, "");
  threeBodTest.calcTBDEnForTor();
  
  int j;
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << threeBodTest.get_energyi_approx(j) << "\t"
    << threeBodTest.get_forcei_approx(j).x();
    cout << ", " << threeBodTest.get_forcei_approx(j).y()<< ", ";
    cout << threeBodTest.get_forcei_approx(j).z()<< "\t"<< endl;
  }
  
  cout << endl;
  
  shared_ptr<ASolver> aSolvall = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                       make_shared<Constants>
                                                       (const_), vals);
  aSolvall->solve_A(1e-4); aSolvall->solve_gradA(1e-4);
  PhysCalc calcphys( aSolvall);
  calcphys.calc_all();
  for ( j = 0; j < num; j++)
  {
    cout << "This is mol " << j << endl;
    cout << calcphys.get_omegai(j) << "\t"
    << calcphys.get_forcei(j)[0];
    cout << ", " << calcphys.get_forcei(j)[1]<< ", ";
    cout << calcphys.get_forcei(j)[2]<< "\t" << endl;
  }
  
  cout << endl;
}

#endif /* ThreeBodyUnitTest_h */
