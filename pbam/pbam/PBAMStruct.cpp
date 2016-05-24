#ifndef __PBAMSTRUCT_H
#define __PBAMSTRUCT_H

//
//  input
//
struct PBAMInput {

  double m_gamma;
  double m_tol;
  double m_pdie;
  double m_sdie;


#ifdef __cplusplus
PBAMInput() :

    m_gamma(0.0001),

  // Solute dielectric
    m_pdie(1.5),

  // Solvent dielectric, from Thomas et. al.
    m_sdie(80),

	   {}
#endif

} ;

//
//  output
//
struct PBAMOutput {

	double m_sumpot,
		m_elecSolvation;

} ;

#endif
