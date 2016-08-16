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
  
  _bessl_consts_ = make_shared<BesselConstants>(2*poles_);
  _bessl_calc_ = make_shared<BesselCalc>(2*poles_, _bessl_consts_);
  _sh_consts_ = make_shared<SHCalcConstants>(2*poles_);
  _sh_calc_ = make_shared<SHCalc>(2*poles_, _sh_consts_);
  _exp_consts_ = make_shared<ExpansionConstants> (poles_);

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
  _setp_ = make_shared<Setup>( 300.0, 0.05, 2., 80., 1, "electrostatics", "tst",
                             false, 100, 0, 15, "tst.map", 1, grid2d,
                             gridax, gridloc, "tst.dx", 1, false, difftype,
                             diffcon, termtype, termval, termnu, confil,
                             conpad, xyzf);
  _syst_ = make_shared<System> ();
  _consts_ = make_shared<Constants> ();
}


PBSAM::PBSAM(string infile)
:
poles_(5),
solveTol_(1e-4)
{
  _setp_ = make_shared<Setup>(infile);
  check_setup();

  _syst_ = make_shared<System> ();
  _consts_ = make_shared<Constants> (*_setp_);

  check_system();
  init_write_system();
  
  _bessl_consts_ = make_shared<BesselConstants>(2*poles_);
  _bessl_calc_ = make_shared<BesselCalc>(2*poles_, _bessl_consts_);
  _sh_consts_ = make_shared<SHCalcConstants>(2*poles_);
  _sh_calc_ = make_shared<SHCalc>(2*poles_, _sh_consts_);
  _exp_consts_ = make_shared<ExpansionConstants> (poles_);
  
  initialize_pbsam();
}


PBSAM::PBSAM(const PBAMInput& pbami, const PBSAMInput& pbsami,
             vector<MoleculeSAM> mls )
:
poles_(5),
solveTol_(1e-4)
{
  int i, j;
  _bessl_consts_ = make_shared<BesselConstants>(2*poles_);
  _bessl_calc_ = make_shared<BesselCalc>(2*poles_, _bessl_consts_);
  _sh_consts_ = make_shared<SHCalcConstants>(2*poles_);
  _sh_calc_ = make_shared<SHCalc>(2*poles_, _sh_consts_);
  _exp_consts_ = make_shared<ExpansionConstants> (poles_);

  // PBSAM part
  vector<string> surffil(pbsami.surfct_),imatfil(pbsami.imatct_);
  vector<string> expfil(pbsami.expct_);
  for (i=0; i<pbsami.surfct_; i++) surffil[i] = string(pbsami.surffil_[i]);
  for (i=0; i<pbsami.imatct_; i++) imatfil[i] = string(pbsami.imatfil_[i]);
  for (i=0; i<pbsami.expct_; i++)   expfil[i] = string(pbsami.expfil_[i]);
  
  // Electrostatics
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

  _setp_ = make_shared<Setup>(pbami.temp_, pbami.salt_, pbami.idiel_,
                             pbami.sdiel_, pbami.nmol_, string(pbami.runType_),
                             string(pbami.runName_), pbami.randOrient_,
                             pbami.boxLen_, pbami.pbcType_, pbami.gridPts_,
                             string(pbami.map3D_), pbami.grid2Dct_,
                             grid2Dname, grid2Dax, grid2Dloc,
                             string(pbami.dxname_), pbami.ntraj_, termcomb,
                             difftype, diffcon, termcond, termval, termnu,
                             confil, conpad, xyzf);
  _setp_->apbs_pbsam_set(surffil, imatfil, expfil);
  check_setup();
  
  vector<shared_ptr<BaseMolecule> > molP(mls.size());
  for (int i = 0; i < mls.size(); i++) molP[i] = make_shared<MoleculeSAM>(mls[i]);
  _syst_ = make_shared<System> (molP, Constants::FORCE_CUTOFF, pbami.boxLen_);
  _consts_ = make_shared<Constants> (*_setp_);
  init_write_system();

  initialize_pbsam();

}

// Function to check inputs from setup file
void PBSAM::check_setup()
{
  try {
    _setp_->check_inputs();
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
    _syst_ = make_shared<System>(*_setp_);
  } catch(const OverlappingMoleculeSAMException& ex1)
  {
    cout << ex1.what() << endl;
    cout << "Provided system has overlapping MoleculeSAMs. ";
    cout << "Please provide a correct system."<< endl;
    exit(0);
  } catch (const NotEnoughCoordsException& ex2)
  {
    cout << ex2.what() << endl;
    exit(0);
  }
  cout << "MoleculeSAM setup okay " << endl;
}


// Rotate MoleculeSAMs if needed and then write out config to pqr
void PBSAM::init_write_system()
{
  if (_setp_->get_randOrient())
  {
    for ( int i = 0; i < _syst_->get_n(); i++)
      _syst_->rotate_mol(i, Quat().chooseRandom());
  }

  // writing initial configuration out
  _syst_->write_to_pqr( _setp_->getRunName() + ".pqr");
  cout << "Written config" << endl;
}

shared_ptr<System> PBSAM::make_subsystem(vector<int> mol_idx)
{
  vector<shared_ptr<BaseMolecule> > sub_mols (mol_idx.size());
  for (int i = 0; i < mol_idx.size(); i++)
  {
    sub_mols[i] = _syst_->get_moli(mol_idx[i]);
  }
  
  shared_ptr<System> _subsys = make_shared<System>(sub_mols,
                                                   _syst_->get_cutoff(),
                                                   _syst_->get_boxlength());
  _subsys -> set_time(_syst_->get_time());
  return _subsys;
}

// Check to see if there are interaction
// matrices provided and if there are expansions from self-polarization
void PBSAM::initialize_pbsam()
{
  int i, k, idx;
  vector<vector<vector<string > > > exp_loc(_syst_->get_n());
  vector<vector<string> > imat_loc(_syst_->get_n());

  cout << "Inside initialize_pbsam" << endl;
  
  for (i = 0; i < _setp_->getNType(); i++)
  {
    string fil=_setp_->getTypeNPQR(i);
    
    if (_setp_->getTypeNImat(i) != "" )
    {
      imat_loc[i].resize(_syst_->get_Ns_i(i));
      string istart = _setp_->getTypeNImat(i);
      imat_loc[i].resize(_syst_->get_Ns_i(i));
      for (int k = 0; k < _syst_->get_Ns_i(i); k++)
        imat_loc[i][k] = istart+to_string(k)+".bin";
    } else
    {
      cout << "Generating IMatrices" << endl;
      // This is done in Solver!
    }
    
    clock_t t3 = clock();

    // Generate surface integrals and buried and exposed points
    // on the molecule surface
    for (k=0; k<_setp_->getTypeNCount(i); k++)
    {
      idx = _syst_->get_mol_global_idx(i,k);
      if (k==0) //Only generate once for each type
      {
        IEMatrix ieMatTest(0, _syst_->get_moli(idx),
                           _sh_calc_, poles_, _exp_consts_, true, 0, true);
        ieMatTest.write_all_mat(fil.substr(0, fil.size()-4));
      } else // copy from first one
      {
        _syst_->copy_grid(_syst_->get_mol_global_idx(i,0), idx);
      }
    }

    
    t3 = clock() - t3;
    printf ("Imat took me %f seconds.\n",
            ((float)t3)/CLOCKS_PER_SEC);
    
    if (_setp_->getTypeNExp(i) != "" )
    {
      string estart = _setp_->getTypeNExp(i);
      exp_loc[i].resize(_syst_->get_Ns_i(i));
      for (int k = 0; k < _syst_->get_Ns_i(i); k++)
      {
        exp_loc[i][k].resize(2);
        exp_loc[i][k][0] = estart+to_string(k) + ".H.exp";
        exp_loc[i][k][1] = estart+to_string(k) + ".F.exp";
      }
    } else if ( _syst_->get_n() != 1 )
    {
      cout << "Solving for self polarization for mol type " << i << endl;
      // Make one MoleculeSAM system for each type to solve for self-polarization
      auto sub_syst = make_subsystem({_syst_->get_mol_global_idx(i,0)});
      Solver self_pol( sub_syst, _consts_, _sh_calc_, _bessl_calc_, poles_);
      self_pol.solve(solveTol_, 100);
      
      //Printing out H and F of selfpol
      //TODO: This should also get printed if system is only 1 Molecule?
      fil = fil.substr(0, fil.size()-4);
      self_pol.get_all_H()[0]->print_all_to_file(_setp_->getTypeNPQR(i)+".H",
                                                 _consts_->get_kappa(),
                                                 _syst_->get_cutoff());
      
      self_pol.get_all_F()[0]->print_all_to_file(fil+".F",
                                                 _consts_->get_kappa(),
                                                 _syst_->get_cutoff());
    }
  }
}


int PBSAM::run()
{
  if ( _setp_->getRunType() == "dynamics")
    run_dynamics();
  else if ( _setp_->getRunType() == "electrostatics")
    run_electrostatics( );
  else if ( _setp_->getRunType() == "energyforce")
    run_energyforce( );
  else if ( _setp_->getRunType() == "bodyapprox")
    run_bodyapprox( );
  else
    cout << "Runtype not recognized! See manual for options" << endl;

  return 0;
}

PBAMOutput PBSAM::run_apbs()
{
  if ( _setp_->getRunType() == "dynamics")
    run_dynamics();
  else if ( _setp_->getRunType() == "electrostatics")
    run_electrostatics( );
  else if ( _setp_->getRunType() == "energyforce")
    run_energyforce( );
  else if ( _setp_->getRunType() == "bodyapprox")
    run_bodyapprox( );
  else
    cout << "Runtype not recognized! See manual for options" << endl;

  PBAMOutput pbsamO;
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
//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, _syst_,
//                                                    _consts_, poles_);
//
//  vector<shared_ptr<BaseTerminate > >  terms(_setp_->get_numterms());
//  for (i = 0; i < _setp_->get_numterms(); i++)
//  {
//    string type = _setp_->get_termtype(i);
//    string bdtype = type.substr(1,1);
//    double val = _setp_->get_termval(i);
//    BoundaryType btype = ( bdtype == "<" ) ? LEQ : GEQ;
//
//    if ( type == "contact" )
//    {
//      cout << "Contact termination found" << endl;
//      double pad = _setp_->get_conpad(j);
//      ContactFile confile (_setp_->get_confile(j));
//      auto conterm = make_shared<ContactTerminate2>(confile, pad);
//
//      terms[i] = make_shared<ContactTerminate2>(confile, pad);
//      j += 1;  // j is index of contact termconditions
//    } else if (type.substr(0,1) == "x")
//    {
//      cout << type << " termination found for MoleculeSAM ";
//      cout << _setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( _setp_->get_termMolIDX(i)[0],
//                                             X, btype, val);
//    } else if (type.substr(0,1) == "y")
//    {
//      cout << type << " termination found for MoleculeSAM ";
//      cout << _setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( _setp_->get_termMolIDX(i)[0],
//                                             Y, btype, val);
//    } else if (type.substr(0,1) == "z")
//    {
//      cout << type << " termination found for MoleculeSAM ";
//      cout << _setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( _setp_->get_termMolIDX(i)[0],
//                                             Z, btype, val);
//    } else if (type.substr(0,1) == "r")
//    {
//      cout << type << " termination found for MoleculeSAM ";
//      cout << _setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
//      terms[i] = make_shared<CoordTerminate>( _setp_->get_termMolIDX(i)[0],
//                                             R, btype, val);
//    } else if (type == "time")
//    {
//      cout << "Time termination found, at time (ps) " << val << endl;
//      terms[i] = make_shared<TimeTerminate>( val);
//    } else cout << "Termination type not recognized!" << endl;
//  }
//
//  cout << "Done making termination conds " << endl;
//  HowTermCombine com = (_setp_->get_andCombine() ? ALL : ONE);
//  auto term_conds = make_shared<CombineTerminate> (terms, com);
//
//  char buff[100], outb[100];
//  sprintf( outb, "%s.stat", _setp_->getRunName().c_str());
//  string statfile = outb;
//
//  for (traj = 0; traj < _setp_->getNTraj(); traj++)
//  {
//    sprintf( buff, "%s_%d.xyz", _setp_->getRunName().c_str(), traj);
//    string xyztraj = buff;
//    sprintf( outb, "%s_%d.dat", _setp_->getRunName().c_str(), traj);
//    string outfile = outb;
//
//    string stats = _setp_->getRunName();
//    _syst_->reset_positions( _setp_->get_trajn_xyz(traj));
//    _syst_->set_time(0.0);
//    BDRun dynamic_run( ASolv, term_conds, outfile);
//    dynamic_run.run(xyztraj, statfile);
//    cout << "Done with trajectory " << traj << endl;
//  }
//
}

void PBSAM::run_electrostatics()
{
  int i;
  clock_t t3 = clock();
  cout << "Before solv constructor" << endl;
  Solver solv( _syst_, _consts_, _sh_calc_, _bessl_calc_, poles_);
  cout << "Before solve" << endl;
  solv.solve(solveTol_, 100);
  
  t3 = clock() - t3;
  printf ("Solve took me %f seconds.\n",
          ((float)t3)/CLOCKS_PER_SEC);
  
  Electrostatic estat(solv.get_all_H(), _syst_, _sh_calc_, _bessl_calc_,
                      _consts_, poles_, _setp_->getGridPts());

  if ( _setp_->getDXoutName() != "" )
    estat.print_dx( _setp_->getDXoutName());

  if ( _setp_->get_3dmap_name() != "" )
    estat.print_3d_heat( _setp_->get_3dmap_name());

  for ( i = 0; i < _setp_->getGridCt(); i++ )
  {
    estat.print_grid(_setp_->getGridAx(i), _setp_->getGridAxLoc(i),
                     _setp_->getGridOutName(i));
  }
}


void PBSAM::run_energyforce()
{
  clock_t t3 = clock();
  Solver solv( _syst_, _consts_, _sh_calc_, _bessl_calc_, poles_);
  solv.solve(solveTol_, 100);
  
  GradSolver gsolv(_syst_, _consts_, _sh_calc_, _bessl_calc_, solv.get_T(),
                   solv.get_all_F(), solv.get_all_H(),
                   solv.get_IE(),
                   solv.get_interpol_list(), _exp_consts_, poles_);
  gsolv.solve(1e-16, 85);
  
  auto focal = make_shared<ForceCalc> (_syst_->get_n(), _syst_->get_all_Ik(),
                                       _consts_->get_dielectric_water(),
                                       _sh_calc_, _bessl_calc_);
  focal->calc_all_f(solv.get_all_H(), solv.get_all_LHN(),
                    gsolv.get_gradH_all(), gsolv.get_gradLHN_all());
  shared_ptr<vector<Pt> > fo = focal->get_all_f();
  
  TorqueCalc tocal(_syst_->get_n());
  tocal.calc_all_tau(_syst_, focal);
  shared_ptr<vector<Pt> > to = tocal.get_all_tau();
  
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
//  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, _syst_,
//                                                    _consts_, poles_);
//  
//  ThreeBody threeBodTest( ASolv, _consts_->get_unitsEnum(), 
//                         _setp_->getRunName(), 100.0);
//  threeBodTest.solveNmer(2);
//  threeBodTest.solveNmer(3);
//  t3 = clock() - t3; 
//  threeBodTest.calcTBDEnForTor();
//
//  threeBodTest.printTBDEnForTor(_setp_->getRunName(), _setp_->getMBDLoc());
//  t3 = clock() - t3;
//  printf ("manybody approx calc took me %f seconds.\n",
//          ((float)t3)/CLOCKS_PER_SEC);
}

