///  @file PBAMFlowWrap.cpp
///  @author  Lisa Felberg
///  @brief c interface for the C++  
///  @ingroup PBAM
///  @version $Id$
///  @attention
///  @verbatim
///
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nathan.baker@pnnl.gov)
///  Pacific Northwest National Laboratory
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright (c) 2010-2015 Battelle Memorial Institute. Developed at the
/// Pacific Northwest National Laboratory, operated by Battelle Memorial
/// Institute, Pacific Northwest Division for the U.S. Department of Energy.
///
/// Portions Copyright (c) 2002-2010, Washington University in St. Louis.
/// Portions Copyright (c) 2002-2010, Nathan A. Baker.
/// Portions Copyright (c) 1999-2002, The Regents of the University of
/// California.
/// Portions Copyright (c) 1995, Michael Holst.
/// All rights reserved.
///
/// Redistribution and use in source and binary forms, with or without
/// modification, are permitted provided that the following conditions are met:
///
/// Redistributions of source code must retain the above copyright notice, this
/// list of conditions and the following disclaimer.
///
/// Redistributions in binary form must reproduce the above copyright notice,
/// this list of conditions and the following disclaimer in the documentation
/// and/or other materials provided with the distribution.
///
/// Neither the name of the developer nor the names of its contributors may be
/// used to endorse or promote products derived from this software without
/// specific prior written permission.
///
/// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
/// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
/// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
/// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
/// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
/// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
/// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
/// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
/// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
/// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
/// THE POSSIBILITY OF SUCH DAMAGE.
///
/// @endverbatim
#include "PBAMWrap.h"
#include "PBAM.h"
#include <iostream>

using namespace std;

struct PBAMInput getPBAMParams()
{
  // create the pbam 
  PBAM pbam;

  // get the struct for use in the c code
  struct PBAMInput pbamI = pbam;
  return pbamI;
}

//
//  for testing only!
//
// struct GeometricFlowOutput runGeometricFlowWrap
//    ( struct GeometricFlowInput geoflowParams )
// {

//    //cout << "boo from GeometricFlowWrap!" << endl; 

//    GeometricFlow GF( geoflowParams );
   
//    AtomList emptyAtomList; // need to fill this with atoms
//    AtomList AL( "imidazole.xyzr", GF.getRadExp(), GF.getFFModel() ); 

//    struct GeometricFlowOutput GFO = GF.run( AL ); //emptyAtomList );
   
//    return GFO;
// }

//
//  print the geometric flow structure for debugging
//
void printPBAMStruct( struct PBAMInput pbamIn )
{
   printf("PBAMInput: %f, %f, %f, %f\n", 
         pbamIn.temp_, 
         pbamIn.idiel_,
         pbamIn.sdiel_,
         pbamIn.salt_);
}

//
//  to call from APBS
//
#ifdef PBAM_APBS
struct PBAMOutput runPBAMWrapAPBS( struct PBAMInput pbamParams,
                                   Valist* molecules[], int nmls ) 
{
  cout << "boo from PBAMWrap!" << endl; 

   // convert Valist to a vector of Molecules
  vector<Molecule> mols;
  cout << "converting atom list" << endl;
  
  //cout << "natoms: " << natoms << endl;
  for (unsigned int mol=0; mol < nmls; mol++) 
  {  
    Vatom *atom;
    unsigned int natoms = Valist_getNumberAtoms(molecules[mol]);

    vector<double> vdw, chg;
    vector<Pt> cgpos;

    for (unsigned int i=0; i < natoms; i++) 
    {   
      atom = Valist_getAtom(molecules[mol], i);
      cout << "i: " << i << endl;

      cgpos.push_back( Pt(Vatom_getPosition(atom)[0],
                          Vatom_getPosition(atom)[1],   
                          Vatom_getPosition(atom)[2])); 
      vdw.push_back(Vatom_getRadius(atom));
      chg.push_back(Vatom_getCharge(atom));
    }
    mols.push_back(Molecule("stat", chg, cgpos, vdw, mol, 0));
  }

  cout << "done with atom list" << endl;
  for (unsigned int mol=0; mol < nmls; mol++) 
  {  
    cout << "This is molecule " << mol << endl;
    for (unsigned int i=0; i < mols[mol].get_m(); i++) 
    { 
      cout << "This is atom " << i << " pos: ";
      cout << mols[mol].get_posj_realspace(i).x() << ", ";
      cout << mols[mol].get_posj_realspace(i).y() << ", ";
      cout << mols[mol].get_posj_realspace(i).z() << endl;
    }
  }  
  
  //  create the PBAM object
  PBAM pbam( pbamParams, mols );

  //  run PBAM!
  struct PBAMOutput pbamO = pbam.run_apbs( );

  return pbamO;
}
#endif // PBAM_APBS