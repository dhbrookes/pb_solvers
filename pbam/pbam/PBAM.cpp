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
poles(5),
solveTol(1e-4)
{
  setp = make_shared<Setup>(infile);
  syst = make_shared<System> ();
  consts = make_shared<Constants> ();
  
  get_check_inputs();

}

// Function to get and check inputs from file
void PBAM::get_check_inputs()
{
  try {
    setp->check_inputs();
  } catch (const BadInputException& ex)
  {
    cout << ex.what() << endl;
    exit(0);
  }
  cout << "All inputs okay " << endl;
  
  consts = make_shared<Constants>(*setp);
  try {
    syst = make_shared<System>(*setp);
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
  
  if (setp->get_randOrient())
  {
    for ( int i = 0; i < syst->get_n(); i++)
      syst->rotate_mol(i, Quat().chooseRandom());
  }
  
  // writing initial configuration out
  syst->write_to_pqr( setp->getRunName() + ".pqr");
  cout << "Written config" << endl;
}



int PBAM::run()
{
  if ( setp->getRunType() == "dynamics")
    run_dynamics();
  else if ( setp->getRunType() == "electrostatics")
    run_electrostatics( );
  else if ( setp->getRunType() == "energyforce")
    run_energyforce( );
  else if ( setp->getRunType() == "bodyapprox")
    run_bodyapprox( );
  else
    cout << "Runtype not recognized! See manual for options" << endl;

  return 0;
}

void PBAM::run_dynamics()
{
  int traj, i, j = 0;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst,
                                                    consts, poles);
  
  vector<shared_ptr<BaseTerminate > >  terms(setp->get_numterms());
  for (i = 0; i < setp->get_numterms(); i++)
  {
    string type = setp->get_termtype(i);
    string bdtype = type.substr(1,2);
    double val = setp->get_termval(i);
    BoundaryType btype = ( bdtype == "<=" ) ? LEQ : GEQ;
    
    if ( type == "contact" )
    {
      cout << "Contact termination found" << endl;
      double pad = setp->get_conpad(j);
      ContactFile confile (setp->get_confile(j));
      auto conterm = make_shared<ContactTerminate2>(confile, pad);

      terms[i] = make_shared<ContactTerminate2>(confile, pad);
      j += 1;  // j is index of contact termconditions
    } else if (type.substr(0,1) == "x")
    {
      cout << type << " termination found for molecule ";
      cout << setp->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp->get_termMolIDX(i)[0],
                                             X, btype, val);
    } else if (type.substr(0,1) == "y")
    {
      cout << type << " termination found for molecule ";
      cout << setp->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp->get_termMolIDX(i)[0],
                                             Y, btype, val);
    } else if (type.substr(0,1) == "z")
    {
      cout << type << " termination found for molecule ";
      cout << setp->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp->get_termMolIDX(i)[0],
                                             Z, btype, val);
    } else if (type.substr(0,1) == "r")
    {
      cout << type << " termination found for molecule ";
      cout << setp->get_termMolIDX(i)[0] << " at a distance " << val << endl;
      terms[i] = make_shared<CoordTerminate>( setp->get_termMolIDX(i)[0],
                                             R, btype, val);
    } else if (type == "time")
    {
      cout << "Time termination found, at time (ps) " << val << endl;
      terms[i] = make_shared<TimeTerminate>( val);
    } else cout << "Termination type not recognized!" << endl;
  }
  
  cout << "Done making termination conds " << endl;
  HowTermCombine com = (setp->get_andCombine() ? ALL : ONE);
  auto term_conds = make_shared<CombineTerminate> (terms, com);
  
  char buff[100], outb[100];
  sprintf( outb, "%s.stat", setp->getRunName().c_str());
  string statfile = outb;
  
  for (traj = 0; traj < setp->getNTraj(); traj++)
  {
    sprintf( buff, "%s_%d.xyz", setp->getRunName().c_str(), traj);
    string xyztraj = buff;
    sprintf( outb, "%s_%d.dat", setp->getRunName().c_str(), traj);
    string outfile = outb;
    
    string stats = setp->getRunName();
    syst->reset_positions( setp->get_trajn_xyz(traj));
    syst->set_time(0.0);
    BDRun dynamic_run( ASolv, term_conds, outfile);
    dynamic_run.run(xyztraj, statfile);
    cout << "Done with trajectory " << traj << endl;
  }


}

void PBAM::run_electrostatics()
{
  int i;
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst,
                                                    consts, poles);
  ASolv->solve_A(solveTol); ASolv->solve_gradA(solveTol);
  Electrostatic Estat( ASolv, setp->getGridPts());
  
  if ( setp->getDXoutName() != "" )
    Estat.print_dx( setp->getDXoutName());
  
  if ( setp->get_3dmap_name() != "" )
    Estat.print_3d_heat( setp->get_3dmap_name());
  
  for ( i = 0; i < setp->getGridCt(); i++ )
  {
    Estat.print_grid(setp->getGridAx(i), setp->getGridAxLoc(i),
                     setp->getGridOutName(i));
  }
}


void PBAM::run_energyforce()
{
  clock_t t;
  t = clock();
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst,
                                                        consts, poles);
  
  ASolv->solve_A(solveTol); ASolv->solve_gradA(solveTol);
  PhysCalc calcEnFoTo( ASolv, setp->getRunName(), consts->get_unitsEnum());
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
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, syst,
                                                    consts, poles);
  ThreeBody threeBodTest( ASolv, consts->get_unitsEnum() );
  threeBodTest.solveNmer(2);
  threeBodTest.solveNmer(3);
  t3 = clock() - t3;
  threeBodTest.calcTBDEnForTor();

  
  threeBodTest.printTBDEnForTor(setp->getMBDLoc());
  
  printf ("manybody approx calc took me %f seconds.\n",
          ((float)t3)/CLOCKS_PER_SEC);

}











