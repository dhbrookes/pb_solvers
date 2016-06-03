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

  vector<double> trans_at_x = {98.33,102.624,94.761,99.071,91.851,83.522,87.029,
    87.205,86.762,93.784,91.864,82.295,78.503,82.203,76.648,80.898,83.62,82.658,
    87.61,92.039,99.156,97.717,98.946,89.004,97.11,89.79,85.491,82.446,88.3,
    88.19,84.849,94.01,88.944,95.098,100.52,101.263};
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
  vector<double> rot_cg_x = {-35.465,-30.388,-29.577,-36.794,-41.493,-30.63,
    -39.732,-49.083,-24.715,-43.218,-40.426,-41.886};
  vector<double> rot_cg_y = {41.47,35.09,40.477,33.798,42.223,41.072,29.147,
    30.885,41.871,33.586,37.415,27.838};
  vector<double> rot_cg_z = {0.47,13.125,-5.638,-6.82,11.224,8.115,9.34,2.646,
    -14.482,7.071,13.57,4.873};

  vector<double> rot1_at_x = {26.345666,26.4538682,21.51591,23.5767556,19.7887928,14.2799782,17.2068573,21.1455619,15.3483584,27.127613,22.8304887,19.0738016,19.6356471,16.847614,8.36684927,11.243107,16.4748788,15.7040622,18.3268588,20.0270005,21.9909676,19.8865605,21.5067207,10.1994736,17.6564658,10.7662936,13.0186869,7.54103477,9.50037956,12.5892817,7.74761095,13.8041871,14.1025636,17.4845116,20.7739425,19.1009159};
  vector<double> rot1_at_y = {55.6408076,58.0642377,48.9047975,49.1148442,45.3031925,39.2734041,41.4981005,38.8779252,45.8247719,43.6317122,48.2636078,39.879489,34.0906375,38.0108013,41.5087489,45.0133433,44.4098445,46.8831797,49.1711484,55.6094971,58.5468531,56.3143925,54.116726,51.246606,50.3099489,47.206651,42.4343654,47.6350506,53.5653759,52.3718688,53.8031408,62.1566371,56.4955772,58.8149954,60.3199327,59.6499494};
  vector<double> rot1_at_z = {16.5345264,8.50589411,4.25129243,-4.00659985,
    2.47901464,4.44861767,4.28204211,6.93151295,9.5228697,8.75867481,14.4709035,
    19.2987094,20.7487153,11.3425847,17.0814019,17.0345497,18.6002767,
    24.6377581,19.6139928,22.2135024,10.5541793,6.20304699,1.66355324,
    2.46292034,-8.15941467,-6.48094454,1.97194765,9.69791467,7.57724451,
    12.0558495,14.7552482,15.8911758,20.8678592,13.873493,
    7.32212478,0.312134665};
  vector<double> rot1_cg_x = {17.2659099,22.6169761,14.6321708,6.84437834,
    22.6556924,23.4835438,12.9687415,6.78620597,11.5642037,13.4379126,
    21.1752768,8.48984003};
  vector<double> rot1_cg_y = {50.9415514,39.1343589,48.2549009,49.3850091,
    52.0267059,44.7913083,43.0756953,52.897983,48.7152127,49.0847765,
    47.4338643,45.2187287};
  vector<double> rot1_cg_z = {9.19753864,16.851665,1.52892834,7.5426126,
    20.2576016,11.5421478,22.1745633,22.929811,-8.1754164,21.3523109,
    22.8018588,20.8884313};

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
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.translate( Pt( 50.0, 50.0, 50.0), 1e48);
  
  for (int i=0; i<molNew.get_nc(); i+=40)
  {
    EXPECT_NEAR( trans_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( trans_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( trans_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<molNew.get_ns(); i+=4)
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
  PQRFile pqr(test_dir_loc + "test_1BRS_cg.pqr");
  Molecule molNew( 0, 0, "stat", pqr.get_charges(),
                  pqr.get_atom_pts(), pqr.get_radii(),
                  pqr.get_cg_centers(), pqr.get_cg_radii());
  molNew.translate( Pt( 100.0, 100.0, 100.0), 90);
  
  for (int i=0; i<molNew.get_nc(); i+=40)
  {
    EXPECT_NEAR( trans_pbc_at_x[ct], molNew.get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( trans_pbc_at_y[ct], molNew.get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( trans_pbc_at_z[ct], molNew.get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<molNew.get_ns(); i+=4)
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
  
  vector<double> sytran_at_x = {9.33366999,13.62767,5.76466999,10.07467,2.85466999,-5.47433001,-1.96733001,-1.79133001,-2.23433001,4.78766999,2.86766999,-6.70133001,-10.49333,-6.79333001,-12.34833,-8.09833001,-5.37633001,-6.33833001,-1.38633001,3.04266999,10.15967,8.72066999,9.94966999,0.00766999287,8.11366999,0.793669993,-3.50533001,-6.55033001,-0.696330007,-0.806330007,-4.14733001,5.01366999,-0.0523300071,6.10166999,11.52367,12.26667};
  vector<double> sytran_at_y = {5.993134,2.507134,-4.921866,-10.738866,-7.977866,-9.068866,-8.605866,-10.039866,-1.637866,-7.564866,0.727133999,-0.752866001,-4.027866,-6.392866,2.521134,3.965134,2.807134,8.637134,6.141134,11.702134,5.635134,2.006134,-2.957866,-0.676866001,-10.594866,-9.382866,-8.039866,2.326134,4.432134,5.439134,9.767134,14.272134,13.427134,9.441134,5.204134,0.816133999};
  vector<double> sytran_at_z = {8.73927798,2.43427798,-0.625722024,-5.10572202,-1.76072202,-1.90972202,-0.894722024,4.33127798,0.113277976,7.84027798,7.51527798,11.251278,14.541278,4.93227798,2.32927798,2.96927798,7.59927798,10.471278,7.90227798,8.65027798,0.842277976,-2.77772202,-4.15972202,-9.87672202,-12.166722,-14.373722,-5.48472202,-5.35672202,-7.54772202,-2.07072202,-3.76872202,-1.88872202,3.63627798,0.175277976,-2.76672202,-8.47472202};
  vector<double> sytran_cg_x = {2.47366999,-3.90633001,1.48066999,-5.19833001,3.22666999,2.07566999,-9.84933001,-8.11133001,2.87466999,-5.41033001,-1.58133001,-11.15833};
  vector<double> sytran_cg_y = {1.065134,-4.011866,-4.822866,2.394134,7.093134,-3.769866,5.332134,14.683134,-9.684866,8.818134,6.026134,7.486134};
  vector<double> sytran_cg_z = {-0.588722024,12.066278,-6.69672202,-7.87872202,10.165278,7.05627798,8.28127798,1.58727798,-15.540722,6.01227798,12.511278,3.81427798};
  
  
  vector<double> sytran1_at_x = {-60.66633,-56.37233,-64.23533,-59.92533,-67.14533,-75.47433,-71.96733,-71.79133,-72.23433,-65.21233,-67.13233,-76.70133,-80.49333,-76.79333,-82.34833,-78.09833,-75.37633,-76.33833,-71.38633,-66.95733,-59.84033,-61.27933,-60.05033,-69.99233,-61.88633,-69.20633,-73.50533,-76.55033,-70.69633,-70.80633,-74.14733,-64.98633,-70.05233,-63.89833,-58.47633,-57.73333};
  vector<double> sytran1_at_y = {-64.006866,-67.492866,-74.921866,-80.738866,-77.977866,-79.068866,-78.605866,-80.039866,-71.637866,-77.564866,-69.272866,-70.752866,-74.027866,-76.392866,-67.478866,-66.034866,-67.192866,-61.362866,-63.858866,-58.297866,-64.364866,-67.993866,-72.957866,-70.676866,-80.594866,-79.382866,-78.039866,-67.673866,-65.567866,-64.560866,-60.232866,-55.727866,-56.572866,-60.558866,-64.795866,-69.183866};
  vector<double> sytran1_at_z = {-61.260722,-67.565722,-70.625722,-75.105722,-71.760722,-71.909722,-70.894722,-65.668722,-69.886722,-62.159722,-62.484722,-58.748722,-55.458722,-65.067722,-67.670722,-67.030722,-62.400722,-59.528722,-62.097722,-61.349722,-69.157722,-72.777722,-74.159722,-79.876722,-82.166722,-84.373722,-75.484722,-75.356722,-77.547722,-72.070722,-73.768722,-71.888722,-66.363722,-69.824722,-72.766722,-78.474722};
  vector<double> sytran1_cg_x = {-67.52633,-73.90633,-68.51933,-75.19833,-66.77333,-67.92433,-79.84933,-78.11133,-67.12533,-75.41033,-71.58133,-81.15833};
  vector<double> sytran1_cg_y = {-68.934866,-74.011866,-74.822866,-67.605866,-62.906866,-73.769866,-64.667866,-55.316866,-79.684866,-61.181866,-63.973866,-62.513866};
  vector<double> sytran1_cg_z = {-70.588722,-57.933722,-76.696722,-77.878722,-59.834722,-62.943722,-61.718722,-68.412722,-85.540722,-63.987722,-57.488722,-66.185722};
  
  vector<double> syrot_at_x = {79.33367,83.62767,75.76467,80.07467,72.85467,64.52567,68.03267,68.20867,67.76567,74.78767,72.86767,63.29867,59.50667,63.20667,57.65167,61.90167,64.62367,63.66167,68.61367,73.04267,80.15967,78.72067,79.94967,70.00767,78.11367,70.79367,66.49467,63.44967,69.30367,69.19367,65.85267,75.01367,69.94767,76.10167,81.52367,82.26667};
  vector<double> syrot_at_y = {75.993134,72.507134,65.078134,59.261134,62.022134,60.931134,61.394134,59.960134,68.362134,62.435134,70.727134,69.247134,65.972134,63.607134,72.521134,73.965134,72.807134,78.637134,76.141134,81.702134,75.635134,72.006134,67.042134,69.323134,59.405134,60.617134,61.960134,72.326134,74.432134,75.439134,79.767134,84.272134,83.427134,79.441134,75.204134,70.816134};
  vector<double> syrot_at_z = {78.739278,72.434278,69.374278,64.894278,68.239278,68.090278,69.105278,74.331278,70.113278,77.840278,77.515278,81.251278,84.541278,74.932278,72.329278,72.969278,77.599278,80.471278,77.902278,78.650278,70.842278,67.222278,65.840278,60.123278,57.833278,55.626278,64.515278,64.643278,62.452278,67.929278,66.231278,68.111278,73.636278,70.175278,67.233278,61.525278};
  vector<double> syrot_cg_x = {72.47367,66.09367,71.48067,64.80167,73.22667,72.07567,60.15067,61.88867,72.87467,64.58967,68.41867,58.84167};
  vector<double> syrot_cg_y = {71.065134,65.988134,65.177134,72.394134,77.093134,66.230134,75.332134,84.683134,60.315134,78.818134,76.026134,77.486134};
  vector<double> syrot_cg_z = {69.411278,82.066278,63.303278,62.121278,80.165278,77.056278,78.281278,71.587278,54.459278,76.012278,82.511278,73.814278};
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
  int nMol = 3; int ct = 0;
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
  
  System sys( mol_, cutoff );
  EXPECT_NEAR( 5.18149787, sys.get_lambda(), preclim);
  EXPECT_NEAR( cutoff/sys.get_cutoff(), 1.0, preclim);

  sys.set_time( 1.435);
  EXPECT_NEAR( 1.435/sys.get_time(), 1.0, preclim);

  EXPECT_EQ(sys.less_than_cutoff(Pt( 5.3, 0.34,  -2.13)), true);
  EXPECT_EQ(sys.less_than_cutoff(Pt(25.3,79.34, -12.13)), false);

  sys.translate_mol(0, Pt(10.0, 3.7, -5.0));
  for (int i=0; i<mol_[0].get_nc(); i+=40)
  {
    EXPECT_NEAR( sytran_at_x[ct], mol_[0].get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( sytran_at_y[ct], mol_[0].get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( sytran_at_z[ct], mol_[0].get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<mol_[0].get_ns(); i+=4)
  {
    EXPECT_NEAR( sytran_cg_x[ct], mol_[0].get_centerk(i).x(), preclim);
    EXPECT_NEAR( sytran_cg_y[ct], mol_[0].get_centerk(i).y(), preclim);
    EXPECT_NEAR( sytran_cg_z[ct], mol_[0].get_centerk(i).z(), preclim);
    ct++;
  }
  
  sys.translate_mol(2, Pt(-10.0, -7.8, -25.0));
  ct = 0;
  for (int i=0; i<mol_[2].get_nc(); i+=40)
  {
    EXPECT_NEAR( sytran1_at_x[ct], mol_[2].get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( sytran1_at_y[ct], mol_[2].get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( sytran1_at_z[ct], mol_[2].get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<mol_[2].get_ns(); i+=4)
  {
    EXPECT_NEAR( sytran1_cg_x[ct], mol_[2].get_centerk(i).x(), preclim);
    EXPECT_NEAR( sytran1_cg_y[ct], mol_[2].get_centerk(i).y(), preclim);
    EXPECT_NEAR( sytran1_cg_z[ct], mol_[2].get_centerk(i).z(), preclim);
    ct++;
  }
  
  sys.rotate_mol(1, Quat( M_PI/2, Pt(0.0, 0.0, 1.0)));
  ct = 0;
  for (int i=0; i<mol_[1].get_nc(); i+=40)
  {
    EXPECT_NEAR( syrot_at_x[ct], mol_[1].get_posj_realspace(i).x(), preclim);
    EXPECT_NEAR( syrot_at_y[ct], mol_[1].get_posj_realspace(i).y(), preclim);
    EXPECT_NEAR( syrot_at_z[ct], mol_[1].get_posj_realspace(i).z(), preclim);
    ct++;
  }
  
  ct = 0;
  for (int i=0; i<mol_[1].get_ns(); i+=4)
  {
    EXPECT_NEAR( syrot_cg_x[ct], mol_[1].get_centerk(i).x(), preclim);
    EXPECT_NEAR( syrot_cg_y[ct], mol_[1].get_centerk(i).y(), preclim);
    EXPECT_NEAR( syrot_cg_z[ct], mol_[1].get_centerk(i).z(), preclim);
    ct++;
  }
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
  double cutoff = 75.0;
  double boxl   = 180.0;
  vector < Molecule > mol_;
  int nMol = 3;
  Pt pos[3] = { Pt(0.0,0.0,0.0), Pt(70.0,70.0,70.0), Pt(-65.3,-68.2,-61.21)};
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
  Pt dis01 = sys.get_pbc_dist_vec_base(sys.get_cogi(0), sys.get_cogi(1));
  EXPECT_NEAR( -70/dis01.x(), 1.0, preclim);
  EXPECT_NEAR( -70/dis01.y(), 1.0, preclim);
  EXPECT_NEAR( -70/dis01.z(), 1.0, preclim);
  
  Pt dis02 = sys.get_pbc_dist_vec_base(sys.get_cogi(0), sys.get_cogi(2));
  EXPECT_NEAR( 65.30/dis02.x(), 1.0, preclim);
  EXPECT_NEAR( 68.20/dis02.y(), 1.0, preclim);
  EXPECT_NEAR( 61.21/dis02.z(), 1.0, preclim);
  
  Pt dis12 = sys.get_pbc_dist_vec_base(sys.get_cogi(1), sys.get_cogi(2));
  EXPECT_NEAR( -44.70/dis12.x(), 1.0, preclim);
  EXPECT_NEAR( -41.80/dis12.y(), 1.0, preclim);
  EXPECT_NEAR( -48.79/dis12.z(), 1.0, preclim);
}


#endif /* SystemUnitTest_h */
