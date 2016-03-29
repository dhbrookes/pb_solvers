//
//  main.cpp
//  pb_solvers_code
//
//  Created by David Brookes on 9/24/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include <memory>
#include "setup.h"
#include "BD.h"
#include "Electrostatics.h"

using namespace std;

int main_dynamics( int poles, double tol, Setup setup,
                  Constants consts, shared_ptr<System> sys)
{
  int i;
  shared_ptr<Constants> constant = make_shared<Constants>(setup);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                    constant, poles);
  
  vector<BaseTerminate> terms(setup.get_numterms());
  for (i = 0; i < setup.get_numterms(); i++)
  {
    string type = setup.get_termtype(i);
    string bdtype = type.substr(1,2);
    double val = setup.get_termval(i);
    BoundaryType btype = ( bdtype == "<=" ) ? LEQ : GEQ;
    
    if ( type == "contact" )
    {
      terms[i] = ContactTerminate( setup.get_termMolIDX(i), val);
    } else if (type.substr(0,1) == "x")
    {
      terms[i] = CoordTerminate( setup.get_termMolIDX(i)[0], X, btype, val);
    } else if (type.substr(0,1) == "y")
    {
      terms[i] = CoordTerminate( setup.get_termMolIDX(i)[0], Y, btype, val);
    } else if (type.substr(0,1) == "z")
    {
      terms[i] = CoordTerminate( setup.get_termMolIDX(i)[0], Z, btype, val);
    } else if (type.substr(0,1) == "r")
    {
      terms[i] = CoordTerminate( setup.get_termMolIDX(i)[0], R, btype, val);
    } else if (type == "time")
    {
      terms[i] = TimeTerminate( val);
    } else cout << "Termination type not recognized!" << endl;
  }
  
  HowTermCombine com = (setup.get_andCombine() ? ALL : ONE);
  
  auto term_conds = make_shared<CombineTerminate> (terms, com);
  BDRun dynamic_run( ASolv, term_conds);
  
  dynamic_run.run();
  
  return 0;
}

int main_electrostatics( int poles, double tol, Setup setup,
                        Constants consts, shared_ptr<System> sys)
{
  int i;
  
  shared_ptr<Constants> constant = make_shared<Constants>(setup);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                    constant, poles);
  ASolv->solve_A(tol); ASolv->solve_gradA(tol);
  Electrostatic Estat( ASolv, setup.getGridPts());
  
  if ( setup.getDXoutName() != "" )
    Estat.print_dx( setup.getDXoutName());
  
  for ( i = 0; i < setup.getGridCt(); i++ )
  {
    Estat.print_grid(setup.getGridAx(i), setup.getGridAxLoc(i),
                     setup.getGridOutName(i));
  }
  return 0;
}


// Main to solve for A and then print energies forces and torques
int main_energyforce( int poles, double tol, Setup setup,
                     Constants consts, shared_ptr<System> sys)
{
  shared_ptr<Constants> constant = make_shared<Constants>(setup);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        constant, poles);
  ASolv->solve_A(tol); ASolv->solve_gradA(tol);
  PhysCalc calcEnFoTo( ASolv, constant->get_unitsEnum());
  calcEnFoTo.calc_all();
  calcEnFoTo.print_all();
  
  return 0;
}

int main_bodyapprox( int poles, double tol, Setup setup,
                  Constants consts, shared_ptr<System> sys)
{
  shared_ptr<Constants> constant = make_shared<Constants>(setup);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                    constant, poles);
  ThreeBody threeBodTest( ASolv, constant->get_unitsEnum() );
  
  threeBodTest.solveNmer(2);
  threeBodTest.solveNmer(3);
  threeBodTest.calcTBDEnForTor();
  threeBodTest.printTBDEnForTor(setup.getMBDLoc());
  
  return 0;
}


int main(int argc, const char * argv[])
{
  string input_file = argv[1];
//  string input_file = "/Users/lfelberg/PBSAM/pb_solvers/pb_solvers_code/test/";
//  input_file += "energyforce_test/run.energyforce.inp";//argv[1];
//  input_file += "electrostatic_test/run.electrostatic.inp";
//  input_file += "dynamics_test/run.dynamics.inp";
//  input_file += "manybodyapprox_test/run.manybodyapprox.inp";

  Setup setp(input_file);
  
  //check inputs:
  try {
    setp.check_inputs();
  } catch (const BadInputException& ex)
  {
    cout << ex.what() << endl;
  }
  cout << "All inputs okay " << endl;

  Constants consts = Constants(setp);
  shared_ptr<System> sys;
  try {
    sys = make_shared<System>(setp);
  } catch(const OverlappingMoleculeException& ex1)
  {
    cout << "Provided system has overlapping molecules. ";
    cout << "Please provide a correct system."<< endl;
  } catch (const NotEnoughCoordsException& ex2)
  {
    cout << ex2.what() << endl;
  }
  cout << "Molecule setup okay " << endl;
  
  if (setp.get_randOrient())
  {
    int i;
    for ( i = 0; i < sys->get_n(); i++)
      sys->rotate_mol(i, Quat().chooseRandom());
  }
  
  // writing initial configuration out
  sys->write_to_pqr( setp.getRunName() + ".pqr");
  cout << "Written config" << endl;
  
  int poles = 10;
  double solv_tol = 1e-4;
  
  if ( setp.getRunType() == "dynamics")
    main_dynamics( poles, solv_tol, setp, consts, sys);
  else if ( setp.getRunType() == "electrostatics")
    main_electrostatics( poles, solv_tol, setp, consts, sys);
  else if ( setp.getRunType() == "energyforce")
    main_energyforce( poles, solv_tol, setp, consts, sys);
  else if ( setp.getRunType() == "bodyapprox")
    main_bodyapprox( poles, solv_tol, setp, consts, sys);
  else
    cout << "Runtype not recognized! See manual for options" << endl;
    
  return 0;
}

