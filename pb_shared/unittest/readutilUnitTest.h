//
//  readutilUnitTest.h
//  pb_solvers_code
//
//  Created by Lisa Felberg on 3/10/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#ifndef readutilUnitTest_h
#define readutilUnitTest_h

#include "readutil.h"

/* Class to test opening PQR and XYZ files */
class ReadUtilUTest : public ::testing::Test
{
public :
  
protected :

  vector<double> sp_x = {2.965, 4.52, 7.092, 7.782, 7.614, 4.82, 4.093, 7.83, 
    6.265, 12.017, 9.993, 8.764, 12.731, 8.46, 4.095, 6.813, 8.422};
  vector<double> sp_z = {-1.084, -3.482, -6.049, 1.145, 4.029, 2.408, -0.923, 
    2.59, -6.542, 3.608, -7.069, 3.456, -4.511, -2.487, -3.01, -0.62, 1.101};
  vector<double> np_y = {-0.636, -0.442, -0.574, -0.594, -0.524, 0.526, -0.048,
   0.512, -0.779, -0.442, 0.201, -0.895, 0.184, -0.595, -0.171, -0.111, -0.75};
  vector<double> np_z = {-0.175, -0.581, -0.81, -0.453, -0.01, 0.656, -0.377,
   0.758, -0.608, 0.723, -0.842, 0.091, -0.82, -0.382, -0.296, -0.774, -0.22};
    
  vector<double> at_x = {48.33,47.142,33.522,41.035,41.864,31.784,30.898,
    33.893,49.156,48.227,39.79,35.538,34.849,43.227,51.263};
  vector<double> at_y = {40.393,25.566,25.331,26.912,35.127,33.217,38.365,
    35.594,40.035,27.504,25.017,34.25,44.167,48.997,35.216};
  vector<double> at_z = {9.798,-1.042,-0.851,4.988,8.574,8.786,4.028,17.666,
    1.901,-6.897,-13.315,-0.595,-2.71,5.886,-7.416};
  vector<double> at_r = {1.824,1.487,1.908,1.387,1.387,1.459,1.824,1.459,1.908,
    1.487,1.908,1.487,1.824,1.487,1.387};
  vector<double> at_c = {0.0966,0.0797,0.5973,0.1007,0.0922,0.1699,-0.4157,
    0.1417,-0.3192,0.0352,-0.0645,0.0791,-0.5163,-0.0122,0.0813};
  
  vector<double> cg_x = {41.47,30.915,44.892,30.335,41.072,27.716,40.051,
    29.456,37.415,23.556};
  vector<double> cg_y = {35.465,33.065,31.202,34.506,30.63,31.582,41.039,
    32.238,40.426,36.176};
  vector<double> cg_z = {0.47,3.732,-3.627,-1.291,8.115,3.988,8.699,15.186,
    13.57,1.911};
  vector<double> cg_r = {13.2578,6.3783,10.5985,4.666,8.7937,4.0027,5.0065,
    3.1349,1.6612,1.6612};
  
  vector<vector<double> > imat = {{},{},{}};
  
  vector<double> scl = {0.21053961,0.21876572, 0.20721523,0.25887287, 0.28581228, 0.22452233, 0.24682216, 0.24901018, 0.23949801, 0.25229589,0.22783195,0.26022692,0.1781166,0.67249496,0.60197448,0.68540096,0.60197448};
  vector<vector<double> > spolHre = {{-0.0569086405,0.0260712491,0.0257008958,0.005235481,-0.0173118192,-0.00856180981,},{-0.00611941969,0.00594485643,0.00287648502,-0.00651224062,-0.00319916935,-0.00294197438,},{0.0437630592,0.00897099853,0.021365849,-0.010795667,0.00812032763,0.00741899161,},{0.0262325702,0.00376868189,-0.00193586639,0.00588471844,-0.00442058501,0.00352281361,},{0.00552877158,0.00343211886,-0.000668543102,0.000670821156,-0.000545964471,-1.62130723e-05,},{0.00389986602,0.00296673714,0.00310364176,0.00140067842,0.00376161814,0.000243639236,},{-0.00775466188,0.00707921297,-0.000750470591,-0.00428955526,-0.00312968678,0.00332338452,},{0.01385662,0.00744574551,0.00599528001,-0.000340534389,0.00271589914,0.00308182702,},{0.00993887881,-0.000675493325,0.00665943185,-0.00595874867,0.00299361641,0.00289549569,},{-0.0149772211,0.00424168113,0.00967218429,0.00315292941,-0.0136369135,-0.0047657203,},{-0.000461199252,0.000431210264,0.00511957525,0.00443801572,-0.00189542418,0.000902068878,},{-0.000249898416,0.0054983672,0.00164742839,-0.00280048924,0.00095732225,0.000312721953,},{-0.011538862,0.00264896378,0.0110225964,0.00602643143,-0.00111149331,-0.00692120288,},{0,0,0,0,0,0,},{-0.00583782653,0.00375597016,0.00211968652,-0.000909173754,-0.00226138201,-8.03077652e-05,},{0,0,0,0,0,0,},{-0.0179288621,0.0124900461,0.00582895032,-0.00504592791,-0.00605247196,-0.0012599674,}};
  vector<vector<double> > spolHim = {{0.0,0.0,0.00964459858,0.0,-0.00757944296,-0.00659596823,},{0.0,0.0,-0.00165224429,0.0,-0.00163459396,0.000665664871,},{0.0,0.0,-0.00799556909,0.0,0.000818580529,-0.00554613199,},{0.0,0.0,-0.00312593082,0.0,-0.00542440823,-0.00804289631,},{0.0,0.0,0.0019811018,0.0,0.00217794045,-3.33506992e-05,},{0.0,0.0,0.00181958663,0.0,0.00245912205,-0.000807222416,},{0.0,0.0,0.00539860387,0.0,-0.00631708316,-0.00593820284,},{0.0,0.0,-0.00378680756,0.0,-0.00194501527,-0.00423593348,},{0.0,0.0,-0.00170593345,0.0,0.00213817344,-0.00135504186,},{0.0,0.0,0.00323553694,0.0,-0.00385667828,-0.00131770137,},{0.0,0.0,0.00110581322,0.0,-0.00128095382,-0.00123417624,},{0.0,0.0,0.000763496041,0.0,0.000432218781,-0.000172166312,},{0.0,0.0,0.00490396708,0.0,0.000829682693,-0.00291276804,},{0.0,0.0,0,0.0,0,0,},{0.0,0.0,0.00202452493,0.0,-0.00215584911,-0.00171813831,},{0.0,0.0,0,0.0,0,0,},{0.0,0.0,0.00159123612,0.0,-0.0018105596,-0.0010116607,}};
  vector<vector<double> > spolFre = {{0.0441857541,-0.0199085448,-0.020878782,-0.00457017258,0.0140744712,0.00798887252,},{0.00366588127,-0.00262835472,-0.000999072376,-4.18060593e-05,0.00127608791,-0.000170762753,},{-0.0718193097,-0.016840448,-0.0381514486,0.0181623461,-0.0143706679,-0.016769726,},{-0.172777903,-0.0563312691,-0.0144647656,-0.0217831081,0.00293096895,-0.00165514847,},{-0.00141046846,-0.000301886988,0.000938720222,0.000600481838,0.000501150525,-0.000572681613,},{-0.00068833022,-0.000602562117,1.47414978e-05,-0.00044595265,2.75120392e-05,5.50230226e-05,},{0.0164583225,-0.012244879,0.00201544505,0.0113395068,0.0119294672,0.00012380472,},{-0.0050765153,-0.00192262841,-0.00257725484,0.000262185395,-0.00115490362,-0.0012712798,},{-0.0489305991,-0.022097321,-0.0102562302,-0.00189178433,-0.00703342908,0.00445122855,},{0.050022961,-0.0348875929,-0.00472978697,0.0172016093,0.00639062461,-0.000670909907,},{-0.0198434657,0.0194169949,-0.00538786785,-0.0126046399,0.00572279708,-0.00518359431,},{0.0145589195,-0.005693701,0.00423613145,-0.000209829042,-0.0016949563,-0.000572450921,},{0.00904922284,0.00294379441,-0.00456172621,-0.00137281415,-0.00289705488,0.00335369683,},{0,0,0,0,0,0,},{0.000685685842,-0.000405964774,-0.000248722109,4.16540907e-05,0.000244034346,-1.93668292e-05,},{0,0,0,0,0,0,},{0.0283733277,-0.0208733848,-0.00907562668,0.0102630238,0.010123643,0.00241235825,}};
  vector<vector<double> > spolFim = {{0.0,0.0,-0.00663657798,0.0,0.00491956015,0.00500337345,},{0.0,0.0,-0.00148021886,0.0,0.00132754129,0.000750911272,},{0.0,0.0,0.0117148143,0.0,0.0020385092,0.0111247663,},{0.0,0.0,-0.0347062376,0.0,-0.0062631668,0.0070352025,},{0.0,0.0,-1.71276303e-06,0.0,-2.61542635e-05,-8.62814537e-05,},{0.0,0.0,-0.0001207205,0.0,-0.000157078511,4.39669436e-05,},{0.0,0.0,-0.00803694108,0.0,0.00630848604,0.0122857756,},{0.0,0.0,0.00137539531,0.0,0.000860872366,0.00151640324,},{0.0,0.0,0.0168670048,0.0,0.00729057417,0.00588229042,},{0.0,0.0,-0.000511724315,0.0,0.00207226658,0.000827838775,},{0.0,0.0,-0.00285877192,0.0,-0.000384615984,-0.000316215871,},{0.0,0.0,0.00538678827,0.0,-0.0029367638,0.00292475064,},{0.0,0.0,-0.000545346304,0.0,-0.000305703867,0.00134106612,},{0.0,0.0,0,0.0,0,0,},{0.0,0.0,-0.000270000994,0.0,0.000265674855,0.000231204042,},{0.0,0.0,0,0.0,0,0,},{0.0,0.0,-0.000646983336,0.0,0.000783653423,0.000237525521,}};
};


// MSMS section
TEST_F(ReadUtilUTest, checkMSMSExceptions)
{
  string path = test_dir_loc + "none.face";
  try
  {
    MSMSFile MSMStest(path);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readMSMS)
{
  string path = test_dir_loc + "test.vert";
  MSMSFile MSMStest(path);
  vector<Pt> my_sp = MSMStest.get_sp();
  vector<Pt> my_np = MSMStest.get_np();

  ASSERT_EQ( my_sp.size(), 1669);
  ASSERT_EQ( my_np.size(), 1669);

  int ct = 0;
  for (int i=0; i<sp_x.size(); i += 100)
  {
    EXPECT_NEAR( my_sp[i].x(), sp_x[ct], preclim);
    EXPECT_NEAR( my_sp[i].z(), sp_z[ct], preclim);
    EXPECT_NEAR( my_np[i].y(), np_y[ct], preclim);
    EXPECT_NEAR( my_np[i].z(), np_z[ct], preclim);
    ct++;
  }

}

// Contact section
TEST_F(ReadUtilUTest, checkContactExceptions)
{
  string path = test_dir_loc + "none.cont";
  try
  {
    ContactFile contTest(path);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readContact)
{
  string path = test_dir_loc + "test.cont";
  ContactFile contTest(path);
  vector<vector<int > > my_pr = contTest.get_at_pairs();

  ASSERT_EQ( contTest.get_moltype1(), 0);
  ASSERT_EQ( contTest.get_moltype2(), 21);

  vector<int > pair1 = {22, 47, 17};
  vector<int > pair2 = {46, 108, 455};
  vector<double > dst = {3.0, 2.57, 9.57};

  for (int i=0; i<my_pr.size(); i++)
  {
    ASSERT_EQ( my_pr[i][0], pair1[i]);
    ASSERT_EQ( my_pr[i][1], pair2[i]);
    EXPECT_NEAR( contTest.get_dists()[i], dst[i], preclim);
  }

}


// TransRot section
TEST_F(ReadUtilUTest, checkTransRotExceptions)
{
  string path = test_dir_loc + "none.transrot";
  try
  {
    PQRFile TRtest(path, 10);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readTransRot)
{
  int M = 3;
  string trrot = test_dir_loc + "test.transrot";
  TransRotFile TRtest(trrot, M);
  vector<double> trx  = { 0.0, 53.20500, -53.20500};
  vector<double> tryy = { 0.0, 92.15376, 92.15376};
  vector<double> trz  = { 0.0, 0.0, 0.0};
  vector<vector<vector<double > > > rot = {
      {{1, 0, 0},              {0, 1, 0},              {0, 0, 1}},
      {{-0.5, -0.866025, 0.0}, { 0.866025, -0.5, 0.0}, {0, 0, 1}},
      {{-0.5,  0.866025, 0.0}, {-0.866025, -0.5, 0.0}, {0, 0, 1}} };
  
  ASSERT_EQ( trrot, TRtest.get_path());
  
  for (int i = 0; i < M; i++)
  {
    EXPECT_NEAR(TRtest.get_trans(i).x(), trx[i], preclim);
    EXPECT_NEAR(TRtest.get_trans(i).y(), tryy[i], preclim);
    EXPECT_NEAR(TRtest.get_trans(i).z(), trz[i], preclim);
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        MyMatrix<double> tmp = TRtest.get_rotmat(i);
        EXPECT_NEAR(tmp(j, k), rot[i][j][k], preclim);
      }
    }
  }
}

// PQR section
TEST_F(ReadUtilUTest, checkPQRExceptions)
{
  string path = test_dir_loc + "none.pqr";
  try
  {
    PQRFile PQRtest(path, 10);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, readPQR)
{
  string PQR = test_dir_loc + "2charge.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  vector<Pt> my_cents = PQRtest.get_cg_centers();
  
  ASSERT_EQ( 1.85, PQRtest.get_cg_radii()[0]);
  ASSERT_EQ( PQR, PQRtest.get_path());
  ASSERT_EQ( my_atoms[0].x(), 0);
  ASSERT_EQ( my_atoms[0].y(), 0);
  ASSERT_EQ( my_atoms[0].z(), 0);
  
  ASSERT_EQ( my_atoms[1].x(), 0);
  ASSERT_EQ( my_atoms[1].y(), 1);
  ASSERT_EQ( my_atoms[1].z(), 0);
  
  ASSERT_EQ(   -1, PQRtest.get_charges()[0]);
  ASSERT_EQ( 1.67, PQRtest.get_charges()[1]);
  ASSERT_EQ( 0.45, PQRtest.get_radii()[0]);
  ASSERT_EQ( 0.20, PQRtest.get_radii()[1]);
  
  ASSERT_EQ(    1, PQRtest.get_cg_centers().size());
  ASSERT_EQ(    0, my_cents[0].x());
  ASSERT_EQ(    0, my_cents[0].y());
  ASSERT_EQ(    0, my_cents[0].z());
  ASSERT_EQ( 1.85, PQRtest.get_cg_radii()[0]);
}

TEST_F(ReadUtilUTest, readPQRNoCen)
{
  string PQR = test_dir_loc + "2charge_nocen.pqr";
  PQRFile PQRtest(PQR, 10);
  vector<Pt> my_atoms = PQRtest.get_atom_pts();
  
  ASSERT_EQ(   4, PQRtest.get_Nc()); 
  ASSERT_EQ( PQR, PQRtest.get_path());
  ASSERT_EQ( my_atoms[0].x(), 0);
  ASSERT_EQ( my_atoms[0].y(), 0);
  ASSERT_EQ( my_atoms[0].z(), 0);
  
  ASSERT_EQ( my_atoms[2].x(), 5);
  ASSERT_EQ( my_atoms[2].y(), 0);
  ASSERT_EQ( my_atoms[2].z(), 0);
  
  ASSERT_EQ(   -1, PQRtest.get_charges()[0]);
  ASSERT_EQ( 7.56, PQRtest.get_charges()[1]);
  ASSERT_EQ( 0.50, PQRtest.get_radii()[0]);
  ASSERT_EQ( 0.87, PQRtest.get_radii()[1]);
}

TEST_F(ReadUtilUTest, readPQRCen)
{
  int ct = 0;
  string PQR = test_dir_loc + "test_1BRS_cg.pqr";
  PQRFile PQRtest(PQR, 2000);
    
  ASSERT_EQ(1403, PQRtest.get_Nc());
  ASSERT_EQ(  47, PQRtest.get_Ns());
  ASSERT_EQ( PQR, PQRtest.get_path());
  
  ct = 0;
  for (int i=0; i < PQRtest.get_Nc(); i+=100)
  {
    EXPECT_NEAR( at_x[ct], PQRtest.get_atom_pts()[i].x(), preclim);
    EXPECT_NEAR( at_y[ct], PQRtest.get_atom_pts()[i].y(), preclim);
    EXPECT_NEAR( at_z[ct], PQRtest.get_atom_pts()[i].z(), preclim);
    EXPECT_NEAR( at_c[ct], PQRtest.get_charges()[i], preclim);
    EXPECT_NEAR( at_r[ct], PQRtest.get_radii()[i], preclim);
    ct++;
  }
}


TEST_F(ReadUtilUTest, checkXYZExceptions)
{
  string path = test_dir_loc + "none.xyz";
  try
  {
    PQRFile XYZtest(path, 10);
    FAIL();
  }
  catch( const CouldNotReadException& err )
  {
    // check exception
    string error_exp = "Could not read: " + path;
    EXPECT_EQ(string(err.what()), error_exp);
  }
}

TEST_F(ReadUtilUTest, checkXYZExceptionShort)
{
  string xyz = test_dir_loc + "2charge_short.xyz";
  try
  {
    XYZFile XYZtest(xyz, 10);
    FAIL();
  }
  catch( const NotEnoughCoordsException& err )
  {
    // check exception
    string error_exp = "File has 5 lines, need 10";
    EXPECT_EQ(string(err.what()), error_exp);
  }
}


TEST_F(ReadUtilUTest, readXYZ)
{
  string xyz = test_dir_loc + "2charge.xyz";
  XYZFile XYZtest(xyz, 10);
  
  ASSERT_EQ(  10, XYZtest.get_nmols());
  ASSERT_EQ( xyz, XYZtest.get_path());
}

TEST_F(ReadUtilUTest, readHFFile)
{
  int p(3), ct;
  
  for( int i = 0; i < 17; i++)
  {
    ct = 0;
    string start = test_dir_loc + "spol_test/test_0.00_p3.0.";
    HFFile htest(start+to_string(i) + ".H.exp", 3);
    HFFile ftest(start+to_string(i) + ".F.exp", 3);
  
    ASSERT_EQ(  3, htest.get_p());
    EXPECT_NEAR(  scl[i], htest.get_kappa(), preclim);
    EXPECT_NEAR(  1000, htest.get_rcut(), preclim);
    ASSERT_EQ( start+to_string(i) + ".H.exp", htest.get_path());
    
    ASSERT_EQ(  3, ftest.get_p());
    EXPECT_NEAR(  scl[i], ftest.get_kappa(), preclim);
    EXPECT_NEAR(  1000, ftest.get_rcut(), preclim);
    ASSERT_EQ( start+to_string(i) + ".F.exp", ftest.get_path());
    
    for (int n = 0; n < p; n++)
    {
      for (int m = 0; m < n+1; m++)
      {
        EXPECT_NEAR( spolHre[i][ct], htest.get_mat_nm(n,m).real(), preclim);
        EXPECT_NEAR( spolHim[i][ct], htest.get_mat_nm(n,m).imag(), preclim);
        EXPECT_NEAR( spolFre[i][ct], ftest.get_mat_nm(n,m).real(), preclim);
        EXPECT_NEAR( spolFim[i][ct], ftest.get_mat_nm(n,m).imag(), preclim);
        ct++;
      }
    }
  }
}


#endif /* readutilUnitTest_h */
