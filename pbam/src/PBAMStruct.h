#ifndef __PBAMSTRUCT_H
#define __PBAMSTRUCT_H

//
//  input
//
struct PBAMInput {

  double temp_;
  double salt_;
  double idiel_;
  double sdiel_;
  char runType_[8192];
  char runName_[8192];

#ifdef __cplusplus
PBAMInput() :
  temp_(298.15),
  salt_(0.01),
  idiel_(1.5), // Solute dielectric
  sdiel_(80.0),
  runType_("energyforce"),
  runName_("tst")
	{}
#endif

} ;

//
//  output
//
struct PBAMOutput {

  double energies_[500];
  double forces_[500][3];

#ifdef __cplusplus
  PBAMOutput() {}
#endif

} ;

#endif
