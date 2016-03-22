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
#include "ThreeBody.h"

using namespace std;

int main_dynamics( int poles, Setup setup, Constants consts, System sys)
{

  
}

int main_electrostatics( int poles, Setup setup, Constants consts, System sys)
{
  
  
}

int main_energyforce( int poles, Setup setup, Constants consts, System sys)
{
  int i, j;
  vector< Molecule > molecules;
  
  for ( i = 0; i < setup.getNType(); i++ )
  
  
  shared_ptr<Constants> const_ = make_shared<Constants>(setup);
  shared_ptr<BesselConstants> bConsta = make_shared<BesselConstants>(2*poles);
  shared_ptr<BesselCalc> bCalcu = make_shared<BesselCalc>(2*poles, bConsta);
  shared_ptr<SHCalcConstants> SHConsta = make_shared<SHCalcConstants>(2*poles);
  shared_ptr<SHCalc> SHCalcu = make_shared<SHCalc>(2*poles, SHConsta);
  shared_ptr<System> sys = make_shared<System>(molecules);
  
  shared_ptr<ASolver> ASolv = make_shared<ASolver> (bCalcu, SHCalcu, sys,
                                                        const_, vals);
  
  PhysCalc calcEnFoTo( ASolv);
  
  
  
}


int main(int argc, const char * argv[])
{
  string input_file = argv[0];
  Setup setup(input_file);

  Constants consts = Constants(setup);
  System sys = System(consts, setup);
  
  int poles = 10;
  
  if ( setup.getRunType() == "dynamics")
    main_dynamics( poles, setup, consts, sys);
  else if ( setup.getRunType() == "potential")
    main_electrostatics( poles, setup, consts, sys);
  else if ( setup.getRunType() == "energyforce")
    main_energyforce( poles, setup, consts, sys)
  else
    cout << "Runtype not recognized! See manual for options" << endl;
    
  return 0;
}






