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
  
  vector<double> cen_map = {9,9,10,21,1,22,22,22,22,20,20,26,35,19,17,31,11,27,
    30,16,9,10,10,8,32,32,1,12,12,2,7,2,7,2,0,10};

  vector<double> trans_at_x = {98.33,102.624,94.761,99.071,91.851,83.522,87.029,87.205,86.762,93.784,91.864,82.295,78.503,82.203,76.648,80.898,83.62,82.658,87.61,92.039,99.156,97.717,98.946,89.004,97.11,89.79,85.491,82.446,88.3,88.19,84.849,94.01,88.944,95.098,100.52,101.263};
  vector<double> trans_at_y = {90.393,86.907,79.478,73.661,76.422,75.331,75.794,74.36,82.762,76.835,85.127,83.647,80.372,78.007,86.921,88.365,87.207,93.037,90.541,96.102,90.035,86.406,81.442,83.723,73.805,75.017,76.36,86.726,88.832,89.839,94.167,98.672,97.827,93.841,89.604,85.216};
  vector<double> trans_at_z = {59.798,53.493,50.433,45.953,49.298,49.149,50.164,55.39,51.172,58.899,58.574,62.31,65.6,55.991,53.388,54.028,58.658,61.53,58.961,59.709,51.901,48.281,46.899,41.182,38.892,36.685,45.574,45.702,43.511,48.988,47.29,49.17,54.695,51.234,48.292,42.584};
  vector<double> trans_cg_x = {91.47,85.09,90.477,83.798,92.223,91.072,79.147,80.885,91.871,83.586,87.415,77.838};
  vector<double> trans_cg_y = {85.465,80.388,79.577,86.794,91.493,80.63,89.732,99.083,74.715,93.218,90.426,91.886};
  vector<double> trans_cg_z = {50.47,63.125,44.362,43.18,61.224,58.115,59.34,52.646,35.518,57.071,63.57,54.873};
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

// TODO: Not sure this is actually something we can test...
TEST_F(MoleculeUTest, checkCreateCen)
{
  srand(1);
  PQRFile pqr(test_dir_loc + "test.pqr");
  MSMSFile surf_file (test_dir_loc + "test.vert");
  vector<Molecule> mols(1);
  mols[0] = Molecule( 0, 0, "stat", pqr.get_charges(),
                      pqr.get_atom_pts(), pqr.get_radii(),
                      surf_file.get_sp(), surf_file.get_np(), 2.5);

//  for (int i=0; i<mols[0].get_ns(); i++)
//  {
//    EXPECT_NEAR( cen_x[i], mols[0].get_centerk(i).x(), preclim);
//    EXPECT_NEAR( cen_y[i], mols[0].get_centerk(i).y(), preclim);
//    EXPECT_NEAR( cen_z[i], mols[0].get_centerk(i).z(), preclim);
//    EXPECT_NEAR( cen_r[i], mols[0].get_ak(i), preclim);
//  }
}

TEST_F(MoleculeUTest, check1BRSCGtoAtMap)
{
  int ct = 0;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  
  vector<Molecule> mols(1);
  mols[0] = Molecule( 0, 0, "stat", pqr.get_charges(),
                      pqr.get_atom_pts(), pqr.get_radii(),
                      pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int i=0; i<mols[0].get_nc(); i+=40)
  {
    EXPECT_EQ( cen_map[ct], mols[0].get_cg_of_ch(i));
    ct++;
  }
}


TEST_F(MoleculeUTest, translate)
{
  int ct = 0;
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.translate( Pt( 50.0, 50.0, 50.0), 1e48);
  
  mols.push_back(molNew);
  System sys(mols);
  sys.write_to_pqr(test_dir_loc+"test_1BRS_trans.pqr");
  
  
  for (int i=0; i<mols[0].get_nc(); i+=40)
  {
    EXPECT_NEAR( trans_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( trans_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( trans_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<mols[0].get_ns(); i+=4)
  {
    EXPECT_NEAR( trans_cg_x[ct], molNew.get_centerk(i).x(), preclim);
    EXPECT_NEAR( trans_cg_y[ct], molNew.get_centerk(i).y(), preclim);
    EXPECT_NEAR( trans_cg_z[ct], molNew.get_centerk(i).z(), preclim);
    ct++;
  }
}

TEST_F(MoleculeUTest, rotateSimple)
{
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  
//  ASSERT_EQ( 0.0, molNew.get_center().x());
//  ASSERT_EQ( 0.0, molNew.get_center().y());
//  ASSERT_EQ( 0.0, molNew.get_center().z());
//  
//  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).x(), preclim);
//  EXPECT_NEAR( 1.0, molNew.get_posj_realspace(0).y(), preclim);
//  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
//  
//  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).x(), preclim);
//  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).y(), preclim);
//  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
//  
//  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
//  
//  EXPECT_NEAR(-1.0, molNew.get_posj_realspace(0).x(), preclim);
//  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).y(), preclim);
//  EXPECT_NEAR( 0.0, molNew.get_posj_realspace(0).z(), preclim);
//  
//  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).x(), preclim);
//  EXPECT_NEAR( -1.0, molNew.get_posj_realspace(1).y(), preclim);
//  EXPECT_NEAR(  0.0, molNew.get_posj_realspace(1).z(), preclim);
}

TEST_F(MoleculeUTest, rotate2)
{
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.rotate( Quat( 1.0, Pt(1.0, 1.0, 1.0)));
  
//  ASSERT_EQ(-10.0, molNew.get_center().x());
//  ASSERT_EQ( 23.4, molNew.get_center().y());
//  ASSERT_EQ( -8.7, molNew.get_center().z());
//  
//  EXPECT_NEAR( 0.693534871, molNew.get_posj(0).x(), preclim);
//  EXPECT_NEAR( 0.639056064, molNew.get_posj(0).y(), preclim);
//  EXPECT_NEAR(-0.332590935, molNew.get_posj(0).z(), preclim);
//  
//  EXPECT_NEAR(-0.332590935, molNew.get_posj(1).x(), preclim);
//  EXPECT_NEAR( 0.693534871, molNew.get_posj(1).y(), preclim);
//  EXPECT_NEAR( 0.639056064, molNew.get_posj(1).z(), preclim);
//  
//  EXPECT_NEAR( -9.30646513, molNew.get_posj_realspace(0).x(), preclim);
//  EXPECT_NEAR( 24.03905610, molNew.get_posj_realspace(0).y(), preclim);
//  EXPECT_NEAR( -9.03259093, molNew.get_posj_realspace(0).z(), preclim);
//  
//  EXPECT_NEAR( -10.3325909, molNew.get_posj_realspace(1).x(), preclim);
//  EXPECT_NEAR( 24.09353490, molNew.get_posj_realspace(1).y(), preclim);
//  EXPECT_NEAR( -8.06094394, molNew.get_posj_realspace(1).z(), preclim);
}

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
