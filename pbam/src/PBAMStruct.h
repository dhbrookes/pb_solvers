#ifndef __PBAMSTRUCT_H
#define __PBAMSTRUCT_H

#define CHR_MAX 1000
#define FIL_MAX 15
#define MOL_MAX 150
#define AT_MAX  50000
#define XYZRCWIDTH 5

//
//  input
//
typedef struct _PBAMInput {

  double temp_;
  double salt_;
  double idiel_;
  double sdiel_;
  int nmol_;
  char runType_[CHR_MAX];
  char runName_[CHR_MAX];
  int randOrient_;
  double boxLen_;
  int pbcType_;


  // Electrostatics
  int gridPts_;
  char map3D_[CHR_MAX];

  int grid2Dct_;
  char grid2D_[FIL_MAX][CHR_MAX];
  char grid2Dax_[FIL_MAX][CHR_MAX];
  double grid2Dloc_[FIL_MAX];

  char dxname_[CHR_MAX];

  // Dynamics
  int ntraj_;
  char termCombine_[CHR_MAX];

  char moveType_[MOL_MAX][CHR_MAX];
  double transDiff_[MOL_MAX];
  double rotDiff_[MOL_MAX];

  int termct_;
  int contct_;

  char termnam_[FIL_MAX][CHR_MAX];
  int termnu_[FIL_MAX][1];
  double termval_[FIL_MAX];
  char confil_[FIL_MAX][CHR_MAX];

  char xyzfil_[MOL_MAX][FIL_MAX][CHR_MAX];
  int xyzct_[MOL_MAX];


 #ifdef __cplusplus
 _PBAMInput() :
   temp_(298.15),
   salt_(0.01),
   idiel_(1.5), // Solute dielectric
   sdiel_(80.0),
   nmol_(1),
   runType_("energyforce"),
   runName_("tst"),
   randOrient_(0),
   boxLen_(1.4e18),
   pbcType_(0),
   gridPts_(15),
   map3D_("tst.map"),
   grid2Dct_(0),
   ntraj_(1),
   termCombine_("or"),
   termct_(1),
   contct_(0)
 	{ }
 #endif

} PBAMInput;

//
//  output
//
typedef struct _PBAMOutput {

  double energies_[MOL_MAX];
  double forces_[MOL_MAX][3];

 #ifdef __cplusplus
   _PBAMOutput() {}
 #endif

} PBAMOutput;

#endif
