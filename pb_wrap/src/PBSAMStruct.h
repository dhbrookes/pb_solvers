#ifndef __PBSAMSTRUCT_H
#define __PBSAMSTRUCT_H

//#define CHR_MAX 1000
//#define FIL_MAX 15
//#define MOL_MAX 150
//#define AT_MAX  50000
//#define XYZRCWIDTH 5

#include "PBAMStruct.h"

//
//  input
//
typedef struct _PBSAMInput {

  char imatfil_[MOL_MAX][CHR_MAX];
  int imatct_;
  
  char surffil_[MOL_MAX][CHR_MAX];
  int surfct_;
  
  char expfil_[MOL_MAX][CHR_MAX];
  int expct_;
  
  double tolsp_;

 #ifdef __cplusplus
 _PBSAMInput() :
  tolsp_(3.0)
 	{ }
 #endif

} PBSAMInput;

//
//  output
//
typedef struct _PBSAMOutput {

  double energies_[MOL_MAX];
  double forces_[MOL_MAX][3];

 #ifdef __cplusplus
   _PBSAMOutput() {}
 #endif

} PBSAMOutput;

#endif
