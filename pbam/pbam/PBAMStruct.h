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

#ifdef __cplusplus
PBAMInput() :

  temp_(298.15),

  salt_(0.01),

  // Solute dielectric
  idiel_(1.5),

  // Solvent dielectric, from Thomas et. al.
  sdiel_(80)

	{}
#endif

} ;

//
//  output
//
struct PBAMOutput {

  double energies[], forces[][3];

} ;

#endif
