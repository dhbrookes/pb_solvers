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

int main_dynamics( Setup setup, Constants consts, System sys)
{
  
  
}

int main_electrostatics( Setup setup, Constants consts, System sys)
{
  
  
}

int main_energyforce( Setup setup, Constants consts, System sys)
{
  
  
}


int main(int argc, const char * argv[])
{
  string input_file = argv[0];
  Setup setup(input_file);

  Constants consts = Constants(setup);
  System sys = System(consts, setup);
  
  if ( setup.getRunType() == "dynamics")
    main_dynamics( consts, sys);
  else if ( setup.getRunType() == "potential")
    main_electrostatics( consts, sys);
  else if ( setup.getRunType() == "energyforce")
    main_energyforce( consts, sys)
  else
    cout << "Runtype not recognized! See manual for options" << endl;
    
  return 0;
}






