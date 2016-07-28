//
//  PBSAM.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 5/24/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "PBSAM.h"

PBSAM::PBSAM() : PBSAMInput()
{
  poles_ = 5;
  solveTol_ = 1e-4;

  vector<string> grid2d = {"tst.2d"};
  vector<string> gridax = {"x"};
  vector<double> gridloc = {0.0};

  // Dynamics
  vector<string> difftype = {"stat"};
  vector<vector <double > > diffcon(1, vector<double> (2));

  vector<string> termtype = {"time"};
  vector<double> termval = {0.0};
  vector<vector<int > > termnu = {{0}};
  vector<string> confil = {""};
  vector<double> conpad = {0.0};
  vector<vector<string > > xyzf = {{""}};


  diffcon[0][0] = 0.0; diffcon[0][1] = 0.0;
  setp_ = make_shared<Setup>( 300.0, 0.05, 2., 80., 1, "electrostatics", "tst",
                             false, 100, 0, 15, "tst.map", 1, grid2d,
                             gridax, gridloc, "tst.dx", 1, false, difftype,
                             diffcon, termtype, termval, termnu, confil,
                             conpad, xyzf);
  syst_ = make_shared<System> ();
  consts_ = make_shared<Constants> ();
}


PBSAM::PBSAM(string infile)
:
poles_(5),
solveTol_(1e-4)
{
  setp_ = make_shared<Setup>(infile);
  check_setup();

  syst_ = make_shared<System> ();
  consts_ = make_shared<Constants> (*setp_);

  check_system();
  init_write_system();
}


PBSAM::PBSAM(const PBSAMInput& pbami, vector<Molecule> mls )
:
poles_(5),
solveTol_(1e-4)
{
  // Electrostatics
  int i, j;
  vector<string> grid2Dname(pbami.grid2Dct_);
  vector<string> grid2Dax(pbami.grid2Dct_);
  vector<double> grid2Dloc(pbami.grid2Dct_);

  for (i=0; i<pbami.grid2Dct_; i++)
  {
    grid2Dname[i] = string(pbami.grid2D_[i]);
    grid2Dax[i] = string(pbami.grid2Dax_[i]);
    grid2Dloc[i] = pbami.grid2Dloc_[i];
  }

  // Dynamics
  bool termcomb = true;
  if (string(pbami.termCombine_) == "or")  termcomb = false;

  vector<string> difftype(pbami.nmol_);
  vector<vector<double > > diffcon(pbami.nmol_, vector<double> (2));
  vector<string> termcond(pbami.termct_);
  vector<double> termval(pbami.termct_);
  vector<vector<int > > termnu(pbami.termct_, vector<int> (1));
  vector<double> conpad;
  vector<string> confil(pbami.contct_);
  vector<vector <string > > xyzf(pbami.nmol_);

  if (string(pbami.runType_) == "dynamics")
  {
    for (i=0; i<pbami.nmol_; i++)
    {
      difftype[i] = string(pbami.moveType_[i]);
      diffcon[i][0] = pbami.transDiff_[i];
      diffcon[i][1] = pbami.rotDiff_[i];
      xyzf[i] = vector<string> (pbami.ntraj_);

      for (j=0; j<pbami.ntraj_; j++)
      {
        xyzf[i][j] = string(pbami.xyzfil_[i][j]);
      }
    }

    for (i=0; i<pbami.termct_; i++)
    {
      termcond[i] = string(pbami.termnam_[i]);
      termval[i] = pbami.termval_[i];
      termnu[i][0] = pbami.termnu_[i][0];

      if (termcond[i] == "contact")
        conpad.push_back(termval[i]);
    }

    for (i=0; i<pbami.contct_; i++)
    {
      confil[i] = string(pbami.confil_[i]);
    }
  }

  setp_ = make_shared<Setup>(pbami.temp_, pbami.salt_, pbami.idiel_,
                             pbami.sdiel_, pbami.nmol_, string(pbami.runType_),
                             string(pbami.runName_), pbami.randOrient_,
                             pbami.boxLen_, pbami.pbcType_, pbami.gridPts_,
                             string(pbami.map3D_), pbami.grid2Dct_,
                             grid2Dname, grid2Dax, grid2Dloc,
                             string(pbami.dxname_), pbami.ntraj_, termcomb,
                             difftype, diffcon, termcond, termval, termnu,
                             confil, conpad, xyzf);

  check_setup();
  
  vector<shared_ptr<Molecule> > molP(mls.size());
  for (int i = 0; i < mls.size(); i++) molP[i] = make_shared<Molecule>(mls[i]);
  
  syst_ = make_shared<System> (molP, Constants::FORCE_CUTOFF, pbami.boxLen_);
  consts_ = make_shared<Constants> (*setp_);
  init_write_system();
}

// Function to check inputs from setup file
void PBSAM::check_setup()
{
  try {
    setp_->check_inputs();
  } catch (const BadInputException& ex)
  {
    cout << ex.what() << endl;
    exit(0);
  }
  cout << "All inputs okay " << endl;
}


// Function to check created system
void PBSAM::check_system()
{
  try {
    syst_ = make_shared<System>(*setp_);
  } catch(const OverlappingMoleculeException& ex1)
  {
    cout << ex1.what() << endl;
    cout << "Provided system has overlapping molecules. ";
    cout << "Please provide a correct system."<< endl;
    exit(0);
  } catch (const NotEnoughCoordsException& ex2)
  {
    cout << ex2.what() << endl;
    exit(0);
  }
  cout << "Molecule setup okay " << endl;
}


// Rotate molecules if needed and then write out config to pqr
void PBSAM::init_write_system()
{
  if (setp_->get_randOrient())
  {
    for ( int i = 0; i < syst_->get_n(); i++)
      syst_->rotate_mol(i, Quat().chooseRandom());
  }

  // writing initial configuration out
  syst_->write_to_pqr( setp_->getRunName() + ".pqr");
  cout << "Written config" << endl;
}



int PBSAM::run()
{
  if ( setp_->getRunType() == "dynamics")
    run_dynamics();
  else if ( setp_->getRunType() == "electrostatics")
    run_electrostatics( );
  else if ( setp_->getRunType() == "energyforce")
    run_energyforce( );
  else if ( setp_->getRunType() == "bodyapprox")
    run_bodyapprox( );
  else
    cout << "Runtype not recognized! See manual for options" << endl;

  return 0;
}

PBSAMOutput PBSAM::run_apbs()
{
  if ( setp_->getRunType() == "dynamics")
    run_dynamics();
  else if ( setp_->getRunType() == "electrostatics")
    run_electrostatics( );
  else if ( setp_->getRunType() == "energyforce")
    run_energyforce( );
  else if ( setp_->getRunType() == "bodyapprox")
    run_bodyapprox( );
  else
    cout << "Runtype not recognized! See manual for options" << endl;

  PBSAMOutput pbsamO;
  return pbsamO;
}

void PBSAM::run_dynamics()
{
//  int traj, i, j = 0;
//  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
//  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
//  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
//  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
//
//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
//                                                    consts_, poles_);
//
//  vector<shared_ptr<BaseTerminate > >  terms(setp_->get_numterms());
//  for (i = 0; i < setp_->get_numterms(); i++)
//  {
//    string type = setp_->get_termtype(i);
//    string bdtype = type.substr(1,1);
//    double val = setp_->get_termval(i);
//    BoundaryType btype = ( bdtype == "<" ) ? LEQ : GEQ;
//
//    if ( type == "contact" )
//    {
//      cout << "Contact termination found" << endl;
//      double pad = setp_->get_conpad(j);
//      ContactFile confile (setp_->get_confile(j));
//      auto conterm = make_shared<ContactTerminate2>(confile, pad);
//
//      terms[i] = make_shared<ContactTerminate2>(confile, pad);
//      j += 1;  // j is index of contact termconditions
//    } else if (type.substr(0,1) == "x")
//    {
//      cout << type << " termination found for molecule ";
//      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
//                                             X, btype, val);
//    } else if (type.substr(0,1) == "y")
//    {
//      cout << type << " termination found for molecule ";
//      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
//                                             Y, btype, val);
//    } else if (type.substr(0,1) == "z")
//    {
//      cout << type << " termination found for molecule ";
//      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
//                                             Z, btype, val);
//    } else if (type.substr(0,1) == "r")
//    {
//      cout << type << " termination found for molecule ";
//      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
//                                             R, btype, val);
//    } else if (type == "time")
//    {
//      cout << "Time termination found, at time (ps) " << val << endl;
//      terms[i] = make_shared<TimeTerminate>( val);
//    } else cout << "Termination type not recognized!" << endl;
//  }
//
//  cout << "Done making termination conds " << endl;
//  HowTermCombine com = (setp_->get_andCombine() ? ALL : ONE);
//  auto term_conds = make_shared<CombineTerminate> (terms, com);
//
//  char buff[100], outb[100];
//  sprintf( outb, "%s.stat", setp_->getRunName().c_str());
//  string statfile = outb;
//
//  for (traj = 0; traj < setp_->getNTraj(); traj++)
//  {
//    sprintf( buff, "%s_%d.xyz", setp_->getRunName().c_str(), traj);
//    string xyztraj = buff;
//    sprintf( outb, "%s_%d.dat", setp_->getRunName().c_str(), traj);
//    string outfile = outb;
//
//    string stats = setp_->getRunName();
//    syst_->reset_positions( setp_->get_trajn_xyz(traj));
//    syst_->set_time(0.0);
//    BDRun dynamic_run( ASolv, term_conds, outfile);
//    dynamic_run.run(xyztraj, statfile);
//    cout << "Done with trajectory " << traj << endl;
//  }
//
}

void PBSAM::run_electrostatics()
{
  int i;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);

//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
//                                                    consts_, poles_);
//  ASolv->solve_A(solveTol_); ASolv->solve_gradA(solveTol_);
//  Electrostatic Estat( ASolv, setp_->getGridPts());
//
//  if ( setp_->getDXoutName() != "" )
//    Estat.print_dx( setp_->getDXoutName());
//
//  if ( setp_->get_3dmap_name() != "" )
//    Estat.print_3d_heat( setp_->get_3dmap_name());
//
//  for ( i = 0; i < setp_->getGridCt(); i++ )
//  {
//    Estat.print_grid(setp_->getGridAx(i), setp_->getGridAxLoc(i),
//                     setp_->getGridOutName(i));
//  }

}


void PBSAM::run_energyforce()
{
  clock_t t3 = clock();
  
//  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
//  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
//  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
//  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
//
//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
//                                                        consts_, poles_);
//
//  ASolv->solve_A(solveTol_); ASolv->solve_gradA(solveTol_);
//  PhysCalc calcEnFoTo( ASolv, setp_->getRunName(), consts_->get_unitsEnum());
//  calcEnFoTo.calc_all();
//  calcEnFoTo.print_all();
  
  t3 = clock() - t3;
  printf ("energyforce calc took me %f seconds.\n",
          ((float)t3)/CLOCKS_PER_SEC);
}

void PBSAM::run_bodyapprox()
{
//  clock_t t3 = clock();  
//  
//  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
//  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
//  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
//  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
//
//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
//                                                    consts_, poles_);
//  
//  ThreeBody threeBodTest( ASolv, consts_->get_unitsEnum(), 
//                         setp_->getRunName(), 100.0);
//  threeBodTest.solveNmer(2);
//  threeBodTest.solveNmer(3);
//  t3 = clock() - t3; 
//  threeBodTest.calcTBDEnForTor();
//
//  threeBodTest.printTBDEnForTor(setp_->getRunName(), setp_->getMBDLoc());
//  t3 = clock() - t3;
//  printf ("manybody approx calc took me %f seconds.\n",
//          ((float)t3)/CLOCKS_PER_SEC);
}

