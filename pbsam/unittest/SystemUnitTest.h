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
  
  vector<double> trans_pbc_at_x = {-31.67,-27.376,-35.239,-30.929,-38.149,-46.478,-42.971,-42.795,-43.238,-36.216,-38.136,-47.705,-51.497,-47.797,-53.352,-49.102,-46.38,-47.342,-42.39,-37.961,-30.844,-32.283,-31.054,-40.996,-32.89,-40.21,-44.509,-47.554,-41.7,-41.81,-45.151,-35.99,-41.056,-34.902,-29.48,-28.737};
  vector<double> trans_pbc_at_y = {50.393,46.907,39.478,33.661,36.422,35.331,35.794,34.36,42.762,36.835,45.127,43.647,40.372,38.007,46.921,48.365,47.207,53.037,50.541,56.102,50.035,46.406,41.442,43.723,33.805,35.017,36.36,46.726,48.832,49.839,54.167,58.672,57.827,53.841,49.604,45.216};
  vector<double> trans_pbc_at_z = {19.798,13.493,10.433,5.953,9.298,9.149,10.164,15.39,11.172,18.899,18.574,22.31,25.6,15.991,13.388,14.028,18.658,21.53,18.961,19.709,11.901,8.281,6.899,1.182,-1.108,-3.315,5.574,5.702,3.511,8.988,7.29,9.17,14.695,11.234,8.292,2.584};
  vector<double> trans_pbc_cg_x = {-38.53,-44.91,-39.523,-46.202,-37.777,-38.928,-50.853,-49.115,-38.129,-46.414,-42.585,-52.162};
  vector<double> trans_pbc_cg_y = {45.465,40.388,39.577,46.794,51.493,40.63,49.732,59.083,34.715,53.218,50.426,51.886};
  vector<double> trans_pbc_cg_z = {10.47,23.125,4.362,3.18,21.224,18.115,19.34,12.646,-4.482,17.071,23.57,14.873};
  
  vector<double> rot_at_x = {-40.393,-36.907,-29.478,-23.661,-26.422,-25.331,-25.794,-24.36,-32.762,-26.835,-35.127,-33.647,-30.372,-28.007,-36.921,-38.365,-37.207,-43.037,-40.541,-46.102,-40.035,-36.406,-31.442,-33.723,-23.805,-25.017,-26.36,-36.726,-38.832,-39.839,-44.167,-48.672,-47.827,-43.841,-39.604,-35.216};
  vector<double> rot_at_y = {48.33,52.624,44.761,49.071,41.851,33.522,37.029,37.205,36.762,43.784,41.864,32.295,28.503,32.203,26.648,30.898,33.62,32.658,37.61,42.039,49.156,47.717,48.946,39.004,47.11,39.79,35.491,32.446,38.3,38.19,34.849,44.01,38.944,45.098,50.52,51.263};
  vector<double> rot_at_z = {9.798,3.493,0.433,-4.047,-0.702,-0.851,0.164,5.39,1.172,8.899,8.574,12.31,15.6,5.991,3.388,4.028,8.658,11.53,8.961,9.709,1.901,-1.719,-3.101,-8.818,-11.108,-13.315,-4.426,-4.298,-6.489,-1.012,-2.71,-0.83,4.695,1.234,-1.708,-7.416};
  vector<double> rot_cg_x = {-35.465,-30.388,-29.577,-36.794,-41.493,-30.63,-39.732,-49.083,-24.715,-43.218,-40.426,-41.886};
  vector<double> rot_cg_y = {41.47,35.09,40.477,33.798,42.223,41.072,29.147,30.885,41.871,33.586,37.415,27.838};
  vector<double> rot_cg_z = {0.47,13.125,-5.638,-6.82,11.224,8.115,9.34,2.646,-14.482,7.071,13.57,4.873};

  vector<double> rot1_at_x = {26.345666,26.4538682,21.51591,23.5767556,19.7887928,14.2799782,17.2068573,21.1455619,15.3483584,27.127613,22.8304887,19.0738016,19.6356471,16.847614,8.36684927,11.243107,16.4748788,15.7040622,18.3268588,20.0270005,21.9909676,19.8865605,21.5067207,10.1994736,17.6564658,10.7662936,13.0186869,7.54103477,9.50037956,12.5892817,7.74761095,13.8041871,14.1025636,17.4845116,20.7739425,19.1009159};
  vector<double> rot1_at_y = {55.6408076,58.0642377,48.9047975,49.1148442,45.3031925,39.2734041,41.4981005,38.8779252,45.8247719,43.6317122,48.2636078,39.879489,34.0906375,38.0108013,41.5087489,45.0133433,44.4098445,46.8831797,49.1711484,55.6094971,58.5468531,56.3143925,54.116726,51.246606,50.3099489,47.206651,42.4343654,47.6350506,53.5653759,52.3718688,53.8031408,62.1566371,56.4955772,58.8149954,60.3199327,59.6499494};
  vector<double> rot1_at_z = {16.5345264,8.50589411,4.25129243,-4.00659985,2.47901464,4.44861767,4.28204211,6.93151295,9.5228697,8.75867481,14.4709035,19.2987094,20.7487153,11.3425847,17.0814019,17.0345497,18.6002767,24.6377581,19.6139928,22.2135024,10.5541793,6.20304699,1.66355324,2.46292034,-8.15941467,-6.48094454,1.97194765,9.69791467,7.57724451,12.0558495,14.7552482,15.8911758,20.8678592,13.873493,7.32212478,0.312134665};
  vector<double> rot1_cg_x = {17.2659099,22.6169761,14.6321708,6.84437834,22.6556924,23.4835438,12.9687415,6.78620597,11.5642037,13.4379126,21.1752768,8.48984003};
  vector<double> rot1_cg_y = {50.9415514,39.1343589,48.2549009,49.3850091,52.0267059,44.7913083,43.0756953,52.897983,48.7152127,49.0847765,47.4338643,45.2187287};
  vector<double> rot1_cg_z = {9.19753864,16.851665,1.52892834,7.5426126,20.2576016,11.5421478,22.1745633,22.929811,-8.1754164,21.3523109,22.8018588,20.8884313};

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


TEST_F(MoleculeUTest, translatePBC)
{
  int ct = 0;
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.translate( Pt( 100.0, 100.0, 100.0), 90);
  
  mols.push_back(molNew);
  System sys(mols);
  sys.write_to_pqr(test_dir_loc+"test_1BRS_trans_pbc.pqr");
  
  
  for (int i=0; i<mols[0].get_nc(); i+=40)
  {
    EXPECT_NEAR( trans_pbc_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( trans_pbc_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( trans_pbc_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<mols[0].get_ns(); i+=4)
  {
    EXPECT_NEAR( trans_pbc_cg_x[ct], molNew.get_centerk(i).x(), preclim);
    EXPECT_NEAR( trans_pbc_cg_y[ct], molNew.get_centerk(i).y(), preclim);
    EXPECT_NEAR( trans_pbc_cg_z[ct], molNew.get_centerk(i).z(), preclim);
    ct++;
  }
}

TEST_F(MoleculeUTest, rotateSimple)
{
  int ct = 0;
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.rotate( Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));

  for (int i=0; i<molNew.get_nc(); i+=40)
  {
    EXPECT_NEAR( rot_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( rot_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( rot_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<molNew.get_ns(); i+=4)
  {
    EXPECT_NEAR( rot_cg_x[ct], molNew.get_centerk(i).x(), preclim);
    EXPECT_NEAR( rot_cg_y[ct], molNew.get_centerk(i).y(), preclim);
    EXPECT_NEAR( rot_cg_z[ct], molNew.get_centerk(i).z(), preclim);
    ct++;
  }
}

TEST_F(MoleculeUTest, rotate2)
{
  int ct = 0;
  vector<Molecule> mols;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.rotate( Quat( 1.0, Pt(1.0, 1.0, 1.0)));
  
  for (int i=0; i<molNew.get_nc(); i+=40)
  {
    EXPECT_NEAR( rot1_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( rot1_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( rot1_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }

  ct = 0;
  for (int i=0; i<molNew.get_ns(); i+=4)
  {
    EXPECT_NEAR( rot1_cg_x[ct], molNew.get_centerk(i).x(), preclim);
    EXPECT_NEAR( rot1_cg_y[ct], molNew.get_centerk(i).y(), preclim);
    EXPECT_NEAR( rot1_cg_z[ct], molNew.get_centerk(i).z(), preclim);
    ct++;
  }
}

class SystemUTest : public ::testing::Test
{
public :
  
protected :
  Constants const_;
  
  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(SystemUTest, checkOverlap)
{
  int M = 2;
  vector < Molecule > mol_;
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int molInd = 0; molInd < M; molInd ++ ) mol_.push_back( molNew );
  mol_[0].translate(Pt(15, 15, 15), 1e48);
  
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
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    mol_.push_back( molNew );
    mol_[molInd].translate(mol_[molInd].get_cog()*(-1), 1e48);
  }
  mol_[0].translate(Pt(15, 15, 15), 1e48);
  mol_[0].translate(Pt(25, 25, 25), 1e48);
  
  try
  {
    System sys( mol_, 80.0, 40.0 );
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
  int nMol = 3;
  double cutoff = 45.876;
  Pt pos[3] = { Pt(0.0,0.0,0.0), Pt(70.0,70.0,70.0), Pt(-70.0,-70.0,-70.0)};
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    mol_.push_back( molNew );
    mol_[molInd].translate(mol_[molInd].get_cog()*(-1)+pos[molInd], 1e48);
  }
  
  // TODO
  System sys( mol_, cutoff );
//  EXPECT_NEAR( 5.7666666667, sys.get_lambda(), preclim);
  EXPECT_NEAR( cutoff/sys.get_cutoff(), 1.0, preclim);

  sys.set_time( 1.435);
  EXPECT_NEAR( 1.435/sys.get_time(), 1.0, preclim);

  EXPECT_EQ(sys.less_than_cutoff(Pt( 5.3, 0.34,  -2.13)), true);
  EXPECT_EQ(sys.less_than_cutoff(Pt(25.3,79.34, -12.13)), false);

  sys.translate_mol(0, Pt(0.0, 0.0, -5.0));
  
//  EXPECT_NEAR(  0.0, sys.get_centeri(0).x(), preclim);
//  EXPECT_NEAR(  0.0, sys.get_centeri(0).y(), preclim);
//  EXPECT_NEAR(-10.0, sys.get_centeri(0).z(), preclim);
//  
//  sys.translate_mol(1, Pt(-10.0, -7.8, -25.0));
//  
//  EXPECT_NEAR(  0.0, sys.get_centeri(1).x(), preclim);
//  EXPECT_NEAR(  0.0, sys.get_centeri(1).y(), preclim);
//  EXPECT_NEAR(  0.0, sys.get_centeri(1).z(), preclim);
//  
//  sys.rotate_mol(1, Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
//  
//  ASSERT_EQ( 0.0, sys.get_centeri(1).x());
//  ASSERT_EQ( 0.0, sys.get_centeri(1).y());
//  ASSERT_EQ( 0.0, sys.get_centeri(1).z());
//  
//  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).x(), preclim);
//  EXPECT_NEAR( 1.0, sys.get_posij(1, 1).y(), preclim);
//  EXPECT_NEAR( 0.0, sys.get_posij(1, 1).z(), preclim);
//  
//  EXPECT_NEAR( -1.0, sys.get_posij(1, 2).x(), preclim);
//  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).y(), preclim);
//  EXPECT_NEAR(  0.0, sys.get_posij(1, 2).z(), preclim);
}

TEST_F(SystemUTest, changeCutoff)
{
  double cutoff = 245.876;
  double boxl   = 255.876;
  vector < Molecule > mol_;
  int nMol = 3;
  Pt pos[3] = { Pt(0.0,0.0,0.0), Pt(70.0,70.0,70.0), Pt(-70.0,-70.0,-70.0)};
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    mol_.push_back( molNew );
    mol_[molInd].translate(mol_[molInd].get_cog()*(-1)+pos[molInd], 1e48);
  }
  
  System sys( mol_, cutoff, boxl );
  EXPECT_NEAR( (boxl/2.0)/sys.get_cutoff(), 1.0, preclim);
}


TEST_F(SystemUTest, PBCcheck)
{
  double cutoff = 125.0;
  double boxl   = 300.0;
  vector < Molecule > mol_;
  int nMol = 3;
  Pt pos[3] = { Pt(0.0,0.0,0.0), Pt(70.0,70.0,70.0), Pt(-70.0,-70.0,-70.0)};
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  
  for (int molInd = 0; molInd < nMol; molInd ++ )
  {
    mol_.push_back( molNew );
    mol_[molInd].translate(mol_[molInd].get_cog()*(-1)+pos[molInd], 1e48);
  }

  System sys( mol_, cutoff, boxl );
  // TODO
  Pt dis01 = sys.get_pbc_dist_vec_base(sys.get_cogi(0), sys.get_cogi(1));
//  EXPECT_NEAR(  10/dis01.x(), 1.0, preclim);
//  EXPECT_NEAR(-7.8/dis01.y(), 1.0, preclim);
//  EXPECT_NEAR(  10/dis01.z(), 1.0, preclim);
  
  Pt dis02 = sys.get_pbc_dist_vec_base(sys.get_cogi(0), sys.get_cogi(2));
//  EXPECT_NEAR( -10/dis02.x(), 1.0, preclim);
//  EXPECT_NEAR( 2.2/dis02.y(), 1.0, preclim);
//  EXPECT_NEAR(   0,     dis02.z(), preclim);
  
  Pt dis12 = sys.get_pbc_dist_vec_base(sys.get_cogi(1), sys.get_cogi(2));
//  EXPECT_NEAR( 0.0,  dis12.x(), preclim);
//  EXPECT_NEAR(  10/dis12.y(), 1.0, preclim);
//  EXPECT_NEAR( -10/dis12.z(), 1.0, preclim);

}


#endif /* SystemUnitTest_h */
