//
//  PBAM.cpp
//  pb_solvers_code
//
//  Created by Lisa Felberg on 5/24/16.
//  Copyright © 2016 Lisa Felberg. All rights reserved.
//

#include "PBAM.h"

PBAM::PBAM(string infile)
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


PBAM::PBAM(const struct PBAMInput& pbami, vector<Molecule> mls )
: 
poles_(5),
solveTol_(1e-4)
{
  setp_ = make_shared<Setup>(pbami.temp_, pbami.salt_, pbami.idiel_,
                             pbami.sdiel_);
  syst_ = make_shared<System> (mls); // TODO: add in boxl and cutoff
  consts_ = make_shared<Constants> (*setp_);
  init_write_system();
}

// Function to check inputs from setup file
void PBAM::check_setup()
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
void PBAM::check_system()
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
void PBAM::init_write_system()
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



int PBAM::run()
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

void PBAM::run_dynamics()
{
  int traj, i, j = 0;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
                                                    consts_, poles_);
  
  vector<shared_ptr<BaseTerminate > >  terms(setp_->get_numterms());
  for (i = 0; i < setp_->get_numterms(); i++)
  {
    string type = setp_->get_termtype(i);
    string bdtype = type.substr(1,2);
    double val = setp_->get_termval(i);
    BoundaryType btype = ( bdtype == "<=" ) ? LEQ : GEQ;
    
    if ( type == "contact" )
    {
      cout << "Contact termination found" << endl;
      double pad = setp_->get_conpad(j);
      ContactFile confile (setp_->get_confile(j));
      auto conterm = make_shared<ContactTerminate2>(confile, pad);

      terms[i] = make_shared<ContactTerminate2>(confile, pad);
      j += 1;  // j is index of contact termconditions
    } else if (type.substr(0,1) == "x")
    {
      cout << type << " termination found for molecule ";
      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
                                             X, btype, val);
    } else if (type.substr(0,1) == "y")
    {
      cout << type << " termination found for molecule ";
      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
                                             Y, btype, val);
    } else if (type.substr(0,1) == "z")
    {
      cout << type << " termination found for molecule ";
      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
                                             Z, btype, val);
    } else if (type.substr(0,1) == "r")
    {
      cout << type << " termination found for molecule ";
      cout << setp_->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp_->get_termMolIDX(i)[0],
                                             R, btype, val);
    } else if (type == "time")
    {
      cout << "Time termination found, at time (ps) " << val << endl;
      terms[i] = make_shared<TimeTerminate>( val);
    } else cout << "Termination type not recognized!" << endl;
  }
  
  cout << "Done making termination conds " << endl;
  HowTermCombine com = (setp_->get_andCombine() ? ALL : ONE);
  auto term_conds = make_shared<CombineTerminate> (terms, com);
  
  char buff[100], outb[100];
  sprintf( outb, "%s.stat", setp_->getRunName().c_str());
  string statfile = outb;
  
  for (traj = 0; traj < setp_->getNTraj(); traj++)
  {
    sprintf( buff, "%s_%d.xyz", setp_->getRunName().c_str(), traj);
    string xyztraj = buff;
    sprintf( outb, "%s_%d.dat", setp_->getRunName().c_str(), traj);
    string outfile = outb;
    
    string stats = setp_->getRunName();
    syst_->reset_positions( setp_->get_trajn_xyz(traj));
    syst_->set_time(0.0);
    BDRun dynamic_run( ASolv, term_conds, outfile);
    dynamic_run.run(xyztraj, statfile);
    cout << "Done with trajectory " << traj << endl;
  }


}

void PBAM::run_electrostatics()
{
  int i;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
                                                    consts_, poles_);
  ASolv->solve_A(solveTol_); ASolv->solve_gradA(solveTol_);
  Electrostatic Estat( ASolv, setp_->getGridPts());
  
  if ( setp_->getDXoutName() != "" )
    Estat.print_dx( setp_->getDXoutName());
  
  if ( setp_->get_3dmap_name() != "" )
    Estat.print_3d_heat( setp_->get_3dmap_name());
  
  for ( i = 0; i < setp_->getGridCt(); i++ )
  {
    Estat.print_grid(setp_->getGridAx(i), setp_->getGridAxLoc(i),
                     setp_->getGridOutName(i));
  }
}


void PBAM::run_energyforce()
{
  clock_t t;
  t = clock();
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
                                                        consts_, poles_);
  
  ASolv->solve_A(solveTol_); ASolv->solve_gradA(solveTol_);
  PhysCalc calcEnFoTo( ASolv, setp_->getRunName(), consts_->get_unitsEnum());
  calcEnFoTo.calc_all();
  calcEnFoTo.print_all();
  t = clock() - t;
  printf ("energyforce calc took me %f seconds.\n",
          ((float)t)/CLOCKS_PER_SEC);

}

void PBAM::run_bodyapprox()
{
  clock_t t3;
  t3 = clock();
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles_);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles_, bConsta);
  auto SHConsta = make_shared<SHCalcConstants>(2*poles_);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles_, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst_,
                                                    consts_, poles_);
  ThreeBody threeBodTest( ASolv, consts_->get_unitsEnum() );
  threeBodTest.solveNmer(2);
  threeBodTest.solveNmer(3);
  t3 = clock() - t3;
  threeBodTest.calcTBDEnForTor();

  
  threeBodTest.printTBDEnForTor(setp_->getMBDLoc());
  
  printf ("manybody approx calc took me %f seconds.\n",
          ((float)t3)/CLOCKS_PER_SEC);

}











