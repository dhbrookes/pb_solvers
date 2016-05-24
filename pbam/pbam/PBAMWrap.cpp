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
  // create the pbam Class
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
   printf("PBAMInput: %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %i, %f\n", 
         pbamIn.m_gamma, 
         pbamIn.m_pdie,
         pbamIn.m_sdie);
}

//
//  to call from APBS
//
#ifdef PBAM_APBS
struct PBAMOutput runPBAMWrapAPBS
   ( struct PBAMInput pbamParams,
     Valist* molecules ) 
{
   //cout << "boo from PBAMWrap!" << endl; 

   //
   //  create the PBAM object
   //
   PBAM pbam( pbamParams );

   //
   // convert Valist to an AtomList
   //
   //cout << "converting atom list" << endl;
   // AtomList atomList;
   // Vatom *atom;
   // unsigned int natoms = Valist_getNumberAtoms(molecules);
   // //cout << "natoms: " << natoms << endl;
   // for (unsigned int i=0; i < natoms; i++) 
   // {		
   //    atom = Valist_getAtom(molecules, i);
   //    //cout << "i: " << i << endl;
   //    Atom myAtom( 
   //          GF.getFFModel(),
   //       Vatom_getPosition(atom)[0],
   //       Vatom_getPosition(atom)[1],		
   //       Vatom_getPosition(atom)[2],		
   //       Vatom_getRadius(atom) * GF.getRadExp(),
   //       Vatom_getCharge(atom) );
   //    atomList.add( myAtom );
   // }
   //cout << "done with atom list" << endl;
   //atomList.print();

   //
   //  run Geoflow!
   //
   struct PBAMOutput pbamO = pbam.run( atomList );
   
   return pbamO;
}
#endif // PBAM_APBS