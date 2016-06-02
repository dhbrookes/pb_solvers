//
//  SystemUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 10/9/15.
//  Copyright © 2015 Lisa Felberg. All rights reserved.
//

#ifndef SystemUnitTest_h
#define SystemUnitTest_h

#include "System.h"

class CGSphereUTest : public ::testing::Test
{
  public :
  
  protected :
  virtual void SetUp() {}
  virtual void TearDown() {}

};

TEST_F(CGSphereUTest, checkUserSpecRadCent)
{
  Pt pos(102.57,33.0,-5.0); double rad = 5.3;
  int M = 3; vector<int> charges(M);
  charges[0]=20; charges[1]=4; charges[2]=1;
  
  CGSphere spNew( pos, rad, charges);
  
  ASSERT_EQ(      3, spNew.get_n());
  ASSERT_EQ(    5.3, spNew.get_a());
  ASSERT_EQ( pos.x(), spNew.get_center().x());
  ASSERT_EQ( pos.y(), spNew.get_center().y());
  ASSERT_EQ( pos.z(), spNew.get_center().z());

  ASSERT_EQ( charges[0], spNew.get_ch()[0]);
  ASSERT_EQ( charges[1], spNew.get_ch()[1]);
  ASSERT_EQ( charges[2], spNew.get_ch()[2]);
}


class MoleculeUTest : public ::testing::Test
{
public :
  
protected :
  virtual void SetUp() {}
  virtual void TearDown() {}

  vector<double> cen_x = {12.1713353,  8.3634097,  6.4158811,  7.4766144,  9.4871234, 11.3504313, 14.7444661, 12.2051735,  9.1604839, 10.0228620, 12.4253020,  9.9452795,  4.1820000,  5.5360000,  5.9370000};
  vector<double> cen_y = {4.7575147,  8.1637544, 12.8213053,  8.2385808,  5.7298088, 11.6602253, 10.8177846, 10.1585197,  1.2516617,  8.6378868,  4.1841707,  1.8740232, 11.6680000,  5.9180000,  3.8990000,};
  vector<double> cen_z = {-0.5078420,  0.7497167, -1.3233676, -3.7726684,  1.7090811,  1.8321553,  4.3286237,  3.4271198,  2.1837470, -4.5740322, -4.1688846,  4.5244564,  0.4430000, -0.7120000,  1.2680000};
  vector<double> cen_r = {4.6373202,  5.0432774,  4.2604095,  4.6023886,  5.2796013,  3.6490073,  2.7461362,  3.8797867,  4.2381621,  3.9138282,  3.7464659,  3.5894635,  1.4590000,  1.6612000,  1.4870000,};

  vector<double> brs_x = {};
  vector<double> brs_y = {};
  vector<double> brs_z = {};
  vector<double> brs_r = {};
};


TEST_F(MoleculeUTest, checkUserSpecCG)
{
  vector<Pt> pos(1); vector<double> cgRad(1);
  int M = 3; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  pos[0] = Pt(2.3, -3.5, 10.1); cgRad[0] = 2.0;
  charges[0]=2.0; vdW[0]=0; posCharges[0] = pos[0];
  charges[1]=2.0; vdW[1]=0; posCharges[1] = pos[0] + Pt(1.0, 0.0, 0.0);
  charges[2]=2.0; vdW[2]=0; posCharges[2] = pos[0] + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( 0, 0, "stat", charges, posCharges, vdW, pos, cgRad, 0, 0);
  
  ASSERT_EQ(      3, molNew.get_nc());
  ASSERT_EQ(      1, molNew.get_ns());
  ASSERT_EQ(    2.0, molNew.get_ak(0));
  ASSERT_EQ( "stat", molNew.get_move_type());
  ASSERT_EQ(    0.0, molNew.get_drot());
  ASSERT_EQ(    0.0, molNew.get_dtrans());

  ASSERT_EQ(      0, molNew.get_cg_of_ch(0));
  ASSERT_EQ(      0, molNew.get_cg_of_ch(1));
  ASSERT_EQ(      0, molNew.get_cg_of_ch(2));
  
  ASSERT_EQ( 2.3, molNew.get_centerk(0).x());
  ASSERT_EQ(-3.5, molNew.get_centerk(0).y());
  ASSERT_EQ(10.1, molNew.get_centerk(0).z());
  
  ASSERT_EQ( 0.0, molNew.get_posj(0).x());
  ASSERT_EQ( 0.0, molNew.get_posj(0).y());
  
  ASSERT_EQ( 1.0, molNew.get_posj(1).x());
  ASSERT_EQ( 0.0, molNew.get_posj(1).z());
  
  ASSERT_EQ( 0.0, molNew.get_posj(2).x());
  ASSERT_EQ( 1.0, molNew.get_posj(2).y());
  
  ASSERT_EQ( 2.3, molNew.get_posj_realspace(0).x());
  ASSERT_EQ(-3.5, molNew.get_posj_realspace(0).y());
  ASSERT_EQ(10.1, molNew.get_posj_realspace(0).z());
  
  ASSERT_EQ( 3.3, molNew.get_posj_realspace(1).x());
  ASSERT_EQ(-3.5, molNew.get_posj_realspace(1).y());
  ASSERT_EQ(10.1, molNew.get_posj_realspace(1).z());
}

TEST_F(MoleculeUTest, checkCreateCen)
{
  srand(1);
  PQRFile pqr(test_dir_loc + "test.pqr");
  MSMSFile surf_file (test_dir_loc + "test.vert");
  vector<Molecule> mols(1);
  mols[0] = Molecule( 0, 0, "stat", pqr.get_charges(),
                      pqr.get_atom_pts(), pqr.get_radii(),
                      surf_file.get_sp(), surf_file.get_np(), 2.5);

  for (int i=0; i<mols[0].get_ns(); i++)
  {
    EXPECT_NEAR( cen_x[i], mols[0].get_centerk(i).x(), preclim);
    EXPECT_NEAR( cen_y[i], mols[0].get_centerk(i).y(), preclim);
    EXPECT_NEAR( cen_z[i], mols[0].get_centerk(i).z(), preclim);
    EXPECT_NEAR( cen_r[i], mols[0].get_ak(i), preclim);
  }
}

TEST_F(MoleculeUTest, checkCreateCen1BRS)
{
  srand(1);
  PQRFile pqr(test_dir_loc + "test_1BRS.pqr");
  MSMSFile surf_file (test_dir_loc + "test_1BRS.vert");
  vector<Molecule> mols(1);
  mols[0] = Molecule( 0, 0, "stat", pqr.get_charges(),
                      pqr.get_atom_pts(), pqr.get_radii(),
                      surf_file.get_sp(), surf_file.get_np(), 3.5);

  EXPECT_NEAR( 2, 2, preclim);
  for (int i=0; i<mols[0].get_ns(); i++)
  {
    // EXPECT_NEAR( brs_x[i], mols[0].get_centerk(i).x(), preclim);
    // EXPECT_NEAR( brs_y[i], mols[0].get_centerk(i).y(), preclim);
    // EXPECT_NEAR( brs_z[i], mols[0].get_centerk(i).z(), preclim);
    // EXPECT_NEAR( brs_r[i], mols[0].get_ak(i), preclim);
  }
}

/*
TEST_F(MoleculeUTest, translate)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, 0, 0, 0.24);
  molNew.translate( Pt( 3.0, -4.5, 10.21), 1e48);
  
  EXPECT_NEAR( -6.5, molNew.get_center().x(), preclim);
  EXPECT_NEAR( 19.4, molNew.get_center().y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_center().z(), preclim);
  
  EXPECT_NEAR( 0.5, molNew.get_posj(0).x(), preclim);
  EXPECT_NEAR(-0.5, molNew.get_posj(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj(0).z(), preclim);
  
  EXPECT_NEAR(-0.5, molNew.get_posj(1).x(), preclim);
  EXPECT_NEAR( 0.5, molNew.get_posj(1).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj(1).z(), preclim);
  
  EXPECT_NEAR( -6.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 18.9, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -7.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( 19.9, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR( 1.51, molNew.get_posj_realspace(1).z(), preclim);
}

TEST_F(MoleculeUTest, rotateSimple)
{
  Pt pos( 0.0, 0.0, 0.0);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, pos, 0, 0, 0.24);
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  ASSERT_EQ( 0.0, molNew.get_center().x());
  ASSERT_EQ( 0.0, molNew.get_center().y());
  ASSERT_EQ( 0.0, molNew.get_center().z());
  
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 1.0, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
  
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  EXPECT_NEAR(-1.0, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
}

TEST_F(MoleculeUTest, rotate2)
{
  Pt pos(-10.0,23.4,-8.7);
  int M = 2; vector<double> charges(M); vector<double> vdW(M);
  vector<Pt> posCharges(M);
  charges[0]=2.0; vdW[0]=3.73; posCharges[0] = pos + Pt(1.0, 0.0, 0.0);
  charges[1]=2.0; vdW[1]=6.32; posCharges[1] = pos + Pt(0.0, 1.0, 0.0);
  
  Molecule molNew( "rot", charges, posCharges, vdW, pos, 0, 0, 0.24);
  molNew.rotate( Quat( 1.0, Pt(1.0, 1.0, 1.0)));
  
  ASSERT_EQ(-10.0, molNew.get_center().x());
  ASSERT_EQ( 23.4, molNew.get_center().y());
  ASSERT_EQ( -8.7, molNew.get_center().z());
  
  EXPECT_NEAR( 0.693534871, molNew.get_posj(0).x(), preclim);
  EXPECT_NEAR( 0.639056064, molNew.get_posj(0).y(), preclim);
  EXPECT_NEAR(-0.332590935, molNew.get_posj(0).z(), preclim);
  
  EXPECT_NEAR(-0.332590935, molNew.get_posj(1).x(), preclim);
  EXPECT_NEAR( 0.693534871, molNew.get_posj(1).y(), preclim);
  EXPECT_NEAR( 0.639056064, molNew.get_posj(1).z(), preclim);
  
  EXPECT_NEAR( -9.30646513, molNew.get_posj_realspace(0).x(), preclim);
  EXPECT_NEAR( 24.03905610, molNew.get_posj_realspace(0).y(), preclim);
  EXPECT_NEAR( -9.03259093, molNew.get_posj_realspace(0).z(), preclim);
  
  EXPECT_NEAR( -10.3325909, molNew.get_posj_realspace(1).x(), preclim);
  EXPECT_NEAR( 24.09353490, molNew.get_posj_realspace(1).y(), preclim);
  EXPECT_NEAR( -8.06094394, molNew.get_posj_realspace(1).z(), preclim);
}
*/
class SystemUTest : public ::testing::Test
{
public :
  
protected :
  Constants const_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
};

/*
TEST_F(SystemUTest, checkOverlap)
{
  vector < Molecule > mol_;
  Pt pos[2] = { Pt( 0.0, 0.0, -5.0), Pt( 0.0, 0.0, 0.0)};
  double rad[2] = { 5.0, 3.7 };
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }
  
  try
  {
    System sys( mol_ );
    FAIL();
  }
  catch( const OverlappingMoleculeException& err )
  {
    // check exception
    string error_exp = "Molecule 0 & 1 overlap";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(SystemUTest, checkPBCOverlap)
{
  vector < Molecule > mol_; const int nMol = 3;
  Pt pos[nMol] = { Pt( 0.0, 0.0, -5.0), Pt(10.0,7.8,25.0), Pt(-10.0,7.8,25.0) };
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0.1; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0.1; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0.1; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", chg, poschg, vdW, pos[molInd], molInd, 0);
    mol_.push_back( molNew );
  }
  
  try
  {
    System sys( mol_, 20.0, 10.0 );
    FAIL();
  }
  catch( const OverlappingMoleculeException& err )
  {
    // check exception
    string error_exp = "Molecule 0 & 1 overlap";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(SystemUTest, checkVals)
{
  vector < Molecule > mol_;
  double cutoff = 45.876;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,7.8,25.0) };
  double rad[3] = { 5.0, 3.7, 8.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }
  
  System sys( mol_, cutoff );
  EXPECT_NEAR( 5.7666666667, sys.get_lambda(), preclim);
  EXPECT_NEAR( cutoff/sys.get_cutoff(), 1.0, preclim);
  
  sys.set_time( 1.435);
  EXPECT_NEAR( 1.435/sys.get_time(), 1.0, preclim);
  
  EXPECT_EQ(sys.less_than_cutoff(Pt( 5.3, 0.34,  -2.13)), true);
  EXPECT_EQ(sys.less_than_cutoff(Pt(25.3,79.34, -12.13)), false);

  sys.translate_mol(0, Pt(0.0, 0.0, -5.0));
  
  EXPECT_NEAR(  0.0, sys.get_centeri(0).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(0).y(), preclim);
  EXPECT_NEAR(-10.0, sys.get_centeri(0).z(), preclim);
  
  sys.translate_mol(1, Pt(-10.0, -7.8, -25.0));
  
  EXPECT_NEAR(  0.0, sys.get_centeri(1).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(1).y(), preclim);
  EXPECT_NEAR(  0.0, sys.get_centeri(1).z(), preclim);
  
  sys.rotate_mol(1, Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
  ASSERT_EQ( 0.0, sys.get_centeri(1).x());
  ASSERT_EQ( 0.0, sys.get_centeri(1).y());
  ASSERT_EQ( 0.0, sys.get_centeri(1).z());
  
  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).x(), preclim);
  EXPECT_NEAR( 1.0, sys.get_posij(1, 1).y(), preclim);
  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).z(), preclim);
  
  EXPECT_NEAR( -1.0, sys.get_posij(1, 2).x(), preclim);
  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).y(), preclim);
  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).z(), preclim);
}

TEST_F(SystemUTest, changeCutoff)
{
  vector < Molecule > mol_;
  double cutoff = 45.876;
  double boxl   = 55.876;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,17.8,15.0) };
  double rad[3] = { 5.0, 3.7, 8.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "stat", rad[molInd], chg, poschg, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }
  
  System sys( mol_, cutoff, boxl );
  EXPECT_NEAR( (boxl/2.0)/sys.get_cutoff(), 1.0, preclim);
}


TEST_F(SystemUTest, PBCcheck)
{
  vector < Molecule > mol_;
  double cutoff =  5.00;
  double boxl   = 20.00;
  Pt pos[3] = { Pt(0.0,0.0,-5.0), Pt(10.0,7.8,25.0), Pt(-10.0,17.8,15.0) };
  double rad[3] = { 1.0, 1.7, 1.6 };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3; vector<double> chg(M); vector<double> vdW(M);
    vector<Pt> poschg(M);
    chg[0]=2.0; vdW[0]=0; poschg[0] = pos[molInd];
    chg[1]=2.0; vdW[1]=0; poschg[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    chg[2]=2.0; vdW[2]=0; poschg[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( "rot", rad[molInd], chg, poschg, vdW, pos[molInd],
                    molInd, 0);
    mol_.push_back( molNew );
  }

  System sys( mol_, cutoff, boxl );
  Pt dis01 = sys.get_pbc_dist_vec(0, 1);
  EXPECT_NEAR(  10/dis01.x(), 1.0, preclim);
  EXPECT_NEAR(-7.8/dis01.y(), 1.0, preclim);
  EXPECT_NEAR(  10/dis01.z(), 1.0, preclim);
  
  Pt dis02 = sys.get_pbc_dist_vec(0, 2);
  EXPECT_NEAR( -10/dis02.x(), 1.0, preclim);
  EXPECT_NEAR( 2.2/dis02.y(), 1.0, preclim);
  EXPECT_NEAR(   0,     dis02.z(), preclim);
  
  Pt dis12 = sys.get_pbc_dist_vec(1, 2);
  EXPECT_NEAR( 0.0,  dis12.x(), preclim);
  EXPECT_NEAR(  10/dis12.y(), 1.0, preclim);
  EXPECT_NEAR( -10/dis12.z(), 1.0, preclim);

}
*/

#endif /* SystemUnitTest_h */
