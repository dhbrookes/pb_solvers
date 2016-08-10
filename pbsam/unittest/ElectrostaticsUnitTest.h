//
//  ElectrostaticsUnitTest.h
//  pbsam
//
//  Created by Felberg, Lisa on 7/27/16.
//  Copyright © 2016 Felberg, Lisa. All rights reserved.
//

#ifndef ElectrostaticsUnitTest_h
#define ElectrostaticsUnitTest_h

#include "Electrostatics.h"

/*
 Class for unit testing energy calculations
 */
class ElectroUTest : public ::testing::Test
{
  protected :
  virtual void SetUp()   {}
  virtual void TearDown() {}
  
  vector<string> dx_head = {"# Data from PBSAM Electrostat run", "# My runname is "+test_dir_loc+"one_mol.dx and units kT", "object 1 class gridpositions counts 10 10 10", "origin -12.1161 -12.8699 -12.9306", "delta 2.37809 0.0e+00 0.0e+00", "delta 0.0e00 2.69144 0.0e+00", "delta 0.0e00 0.0e+00 2.35765", "object 2 class gridconnections counts 10 10 10", "object 3 class array type double rank 0 items 1000 data follows"};
  
  vector<double> dx_vals = {-0.016893876,-0.018298784,-0.019045288,-0.018782092,-0.017286529,-0.014621549,-0.011185659,-0.007582257,-0.004371536,-0.001873302};
  
  vector<vector<double> > map3d_3 = {{-2.736622041,-0.221695402,-7.262395704,-0.000558689},{-0.882927438,-8.891490401,-0.048928185,0.000011763},{6.482287959,-1.807057062,2.651447472,0.000210335},{8.219685149,3.310386044,3.606237079,-0.000003777},{-6.386156024,4.704504578,2.289000672,-0.000066206},{-5.159104456,4.966351785,-5.147400762,-0.000279652},{-6.956681763,3.306713797,-0.054286748,-0.000172178},{2.643386044,0.708623303,5.870850177,0.000112214},{-0.526065784,-3.338715944,6.654473614,-0.000007009},{4.528601076,1.457047585,-7.377275489,-0.000017660},{5.875638168,-5.439137567,-6.685115631,0.000019921},{6.677545244,-0.713647957,-2.252582765,0.000068430},{4.150879665,1.969677590,-2.739859042,-0.000053576},{-4.074101286,0.072572962,-6.927044006,-0.000580607},{8.416373064,1.585883378,15.526008411,-0.000085804},{18.472032799,0.397811147,18.228702188,0.000204220},{15.186729729,14.499471439,18.856579361,0.000080592},{6.970098031,9.987671027,16.534792168,-0.000084379},{9.758847742,15.081184965,10.425959445,-0.000093421},{5.165443086,15.436882280,15.519073459,-0.000013128},{17.920100998,6.984521697,20.972127869,0.000191270},{17.018670767,1.873126262,22.227929855,0.000112044},{11.670181450,6.632060483,8.078605783,-0.000480499},{14.566561389,1.600499565,8.888757235,-0.000051098},{20.381967111,1.392157210,13.321290102,0.000154262},{13.165378810,14.850734893,17.156710819,0.000049990},{10.272244749,9.597925960,8.717774486,-0.000466293},{-16.040261613,-15.677019212,-16.732599674,-0.000082631},{-7.410059408,-13.497417018,-14.077442472,0.000195238},{-7.246120290,0.378421995,-12.819770421,-0.000101103},{-9.985258627,-3.609801487,-8.940075204,0.000055998},{-18.270907809,-0.704331254,-21.280120313,-0.000199049},{-17.900041355,3.443602743,-16.157566333,-0.000065474},{-9.620963799,-8.591393229,-11.117619683,0.000185732},{-14.318594888,-15.057048079,-9.448613903,-0.000025187},{-15.266605152,-3.248900689,-23.597865454,-0.000250762},{-8.232721105,-8.821783535,-22.934142510,-0.000041543},{-11.906610375,-14.948505843,-18.354837032,-0.000053533},{-12.548602030,-0.407878134,-14.873427283,-0.000000328},{-17.499906355,-5.810827923,-23.062709534,-0.000385286}};
  
};

TEST_F(ElectroUTest, one_mol_test)
{
  int pol(3), nmol(1);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<MoleculeSAM> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> (kT);
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_MoleculeSAM(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "spol_test/test_0.00_p3.0.";
  vector<vector<string> > imat_loc(sys->get_n());
  string exp_str;
  
  vector<shared_ptr<HMatrix> > Hmat(sys->get_n());
  
  for (int i = 0; i < sys->get_n(); i++)
  {
    Hmat[i] = make_shared<HMatrix>(i, sys->get_Ns_i(i), pol, cst->get_kappa());
    imat_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      imat_loc[i][k] = istart+to_string(k) + ".out.bin";
      exp_str = estart+to_string(k)+ ".H.exp";
      Hmat[i]->init_from_exp(exp_str, k);
    }
  }

  string dxmap = "one_mol.dx";
  Electrostatic estat(Hmat, sys, SHCalcTest, BesselCal, cst, pol, 10);
  estat.print_dx(test_dir_loc + "one_mol.dx");
  
  string inputLine;
  ifstream fin(test_dir_loc+dxmap);
  getline(fin,inputLine);
  
  int ct(0), offset(9);
  while ((ct < offset) && (!fin.eof()))
  {
    EXPECT_TRUE(inputLine == dx_head[ct]);
    getline(fin,inputLine);
    ct++;
  }
  
  int inn(0);
  ct = 0;
  vector<vector<vector<double > > > potdx = estat.get_potential();
  for(int i = 0; i<estat.get_bins().size(); i++)
    for(int j = 0; j<estat.get_bins().size(); j++)
      for(int k = 0; k<estat.get_bins().size(); k++)
      {
        if (ct % 100 == 0)
        {
          EXPECT_NEAR(potdx[i][j][k], dx_vals[inn], preclim);
          inn++;
        }
        ct++;
      }
}

TEST_F(ElectroUTest, three_mol_test)
{
  int pol(3), nmol(3);
  PQRFile pqr(test_dir_loc + "test_cged.pqr");
  vector<shared_ptr<MoleculeSAM> > mols;
  for (int i=0; i<nmol; i++)
    mols.push_back(make_shared<MoleculeSAM>(0, 0, "stat", pqr.get_charges(),
                                         pqr.get_atom_pts(), pqr.get_radii(),
                                         pqr.get_cg_centers(),
                                         pqr.get_cg_radii()));
  
  mols[0]->translate(Pt(-9.28786458,-7.35779167,-0.15628125), 1e14);
  mols[1]->translate(Pt(3.71213542,-0.35779167,14.84371875), 1e14);
  mols[2]->translate(Pt(-22.28786458,-14.35779167,-15.15628125), 1e14);
  auto sys = make_shared<System>(mols);
  auto cst = make_shared<Constants> ();
  cst->set_dielectric_water(80);
  cst->set_dielectric_prot(4);
  cst->set_salt_concentration(0.01);
  cst->set_temp(298.15);
  cst->set_kappa(0.0325628352);
  
  auto _SHConstTest = make_shared<SHCalcConstants> (2*pol);
  auto SHCalcTest = make_shared<SHCalc> (2*pol, _SHConstTest);
  auto BesselCons = make_shared<BesselConstants> (2*pol);
  auto BesselCal = make_shared<BesselCalc>(2*pol, BesselCons);
  auto _expcons = make_shared<ExpansionConstants> (pol);
  
  // Generate surface integrals
  for (int i=0; i<nmol; i++)
    IEMatrix ieMatTest(0, sys->get_MoleculeSAM(i),
                       SHCalcTest, pol, _expcons, true, 0, true);
  
  string istart = test_dir_loc + "imat_test/imat.sp";
  string estart = test_dir_loc + "grad_test/mpol3.";
  vector<vector<string> > imat_loc(sys->get_n());
  vector<vector<string > > exp_loc(sys->get_n());
  vector<shared_ptr<HMatrix> > Hmat(sys->get_n());
  
  // Generate surface integrals
  for (int i = 0; i < sys->get_n(); i++)
  {
    Hmat[i] = make_shared<HMatrix>(i, sys->get_Ns_i(i), pol, cst->get_kappa());
    imat_loc[i].resize(sys->get_Ns_i(i));
    exp_loc[i].resize(sys->get_Ns_i(i));
    for (int k = 0; k < sys->get_Ns_i(i); k++)
    {
      imat_loc[i][k] = istart+to_string(k)+ ".out.bin";
      exp_loc[i][k] = estart+to_string(i)+"."+to_string(k)+".H.exp";
      Hmat[i]->init_from_exp(exp_loc[i][k], k);
    }
  }
  
  string map3d = "threemol_map.out";
  Electrostatic estat(Hmat, sys, SHCalcTest, BesselCal, cst, pol, 2);
  estat.print_3d_heat(test_dir_loc + map3d);
  
  string inputLine;
  ifstream fin(test_dir_loc+map3d);
  getline(fin,inputLine);
  double px, py, pz, pot;
  
  int ct(0), inn(0), offset(5);
  while ((!fin.eof()))
  {
    if ((ct >= offset) && ((ct-offset)%100 == 0))
    {
      istringstream ss( inputLine );
      ss >> px >> py >> pz >> pot;
      EXPECT_NEAR( px, map3d_3[inn][0], preclim);
      EXPECT_NEAR( py, map3d_3[inn][1], preclim);
      EXPECT_NEAR( pz, map3d_3[inn][2], preclim);
      EXPECT_NEAR( pot, map3d_3[inn][3], preclim);
      inn++;
    }
    getline(fin,inputLine);
    ct++;
  }
}


#endif /* ElectrostaticsUnitTest_h */
