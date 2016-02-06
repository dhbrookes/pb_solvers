//
//  ASolverUnitTest.h
//  pbsolvers
//
//  Created by Marielle Soniat on 10/6/15.
//  Copyright (c) 2015 Marielle Soniat. All rights reserved.
//

#ifndef pbsolvers_ASolverUnitTest_h
#define pbsolvers_ASolverUnitTest_h

#include "ASolver.h"

class ASolverUTest : public ::testing::Test
{
public :
  
protected :
  
  int vals_;
  Constants const_;
  vector< Molecule > mol_;
  vector< Molecule > mol_sing_;
  
  virtual void SetUp()
  {
    mol_.clear( );
    Pt pos[2]     = { Pt( 0.0, 0.0, -5.0 ), Pt( 10.0, 7.8, 25.0 ) };
    Pt cgPos[2]   = { Pt( 0.0, 0.0, -5.5 ), Pt( 11.0, 6.9, 24.3 ) };
    double cg[2] = { 5.0, -0.4};
    double rd[2] = { 5.6, 10.4};
    
    Pt cgPosSi[2] = { Pt( 0.0, 0.0, -35.0 ), Pt( 0.0, 0.0, 0.0 ) };
    
    for (int molInd = 0; molInd < 2; molInd ++ )
    {
      int M = 1;
      vector<double> charges(1);
      vector<Pt> posCharges(1);
      
      charges[0]    = cg[molInd];
      posCharges[0] = cgPos[molInd];
      
      Molecule molNew( M, rd[molInd], charges, posCharges, pos[molInd]);
      mol_.push_back( molNew );
      
      charges[0]    = 2.0;
      posCharges[0] = cgPosSi[molInd];

      Molecule molSing( M, 10.0, charges, posCharges);
      mol_sing_.push_back( molSing );
    }
  } // end SetUp
  
  virtual void TearDown() {}
  
  double ATrip0[15]   = {6.01094069,0.0141280225,1.03433512,-0.804854118,
    1.24302352e-05,2.58107532e-05,5.8366607e-05,-0.182460812,3.82002752e-06,
    0.23555681,0.16212821,-7.82709e-09,3.914100e-07,2.5549544e-08,0.22607375};
  double ATrip0im[15] = {0,0,1.03698189,0,0.000412418784,3.23186115e-06,0,
    -0.18242449,3.11917432e-07,-0.235553474,0,2.63554468e-06,1.04180458e-08,
    5.16522503e-07,1.02665424e-08};
  double ATrip1[15]   = {6.01119603,-0.00735850845,1.01433618,-0.80724428,
    0.000252370073,0.00231355778,-3.33699295e-05,-0.182293154,-2.53839549e-06,
    0.23531127,0.16214575,1.7126064e-06,-1.9023086e-05,-9.1676707e-08,
    0.22609914};
  double ATrip1im[15] = {0,0,1.03337432,0,0.000193710438,-4.13320113e-05,0,
    -0.182483415,-9.28886031e-06,-0.235542082,0,1.32479048e-06,1.79437885e-06,
    2.38805876e-07,-2.11951326e-06};
  double ATrip2[15]   = {6.01116143,-0.00712987242,1.05234476,-0.807013253,
    -0.000264824695,0.0019962685,-3.00198469e-05,-0.182594503,-3.60529882e-06,
    0.235757586,0.162141472,-1.7044803e-06,-1.46864088e-05,6.58681972e-08,
    0.22609354};
  double ATrip2im[15] = {0,0,1.03341635,0,0.000184134752,3.80938095e-05,0,
    -0.182482137,9.60094345e-06,-0.235542368,0,1.19249837e-06,-1.80475738e-06,
    2.8073034e-07,2.10915172e-06};

  double A0[15] = {5.06502332,-0.459077988,-0.000488586514,0.0313638895,
    -0.0001113252,-2.2608511e-05,-0.00206626713,-1.20922746e-05,-5.8562079e-06,
    -9.08691077e-07,0.000132860057,-3.86433215e-07,-9.97731223e-07,
    -3.82176443e-07,-3.75510255e-08};
  double A0_im[15] = {0,0,-0.00021783804,0,-3.85223997e-05,-2.08028337e-05,0,
    -8.64506194e-07,-2.89571751e-06,-1.2658108e-06,0,9.07199987e-07,
    2.00374768e-07,-2.07633434e-07,-8.16687899e-08};
  double A1[15] = {-0.396665696,-0.157791465,-0.100108147,0.113721012,
    0.0563843732,0.00185962096,-0.0415378387,-0.0286835422,-0.00276725208,
    0.00189910059,0.0110646934,0.013544291,0.00205846384,-0.00116220012,
    -5.31391954e-06};
  double A1_im[15] = {0,0,0.00922544874,0,0.0293832761,0.0232705597,0,
    -0.0222273022,-0.0142638165,-0.000444043694,0,0.0109506821,0.00839150215,
    0.00171418682,0.000215801105};
  
  double AMul0Sing[15]   = {6.0131921,0.451572318,1.06409876,-0.565532819,
    0.0265785191,0,0.107782872,-0.166377441,0,0.235821288,0.208159863,
    0.00831875258,0,0.000298201925,0.226117001};
  double AMul0Singim[15] = {0,0,1.06409876,0,0.0265785191,0,0,-0.166377441,0,
    -0.235821288,0,0.00831875258,0,-0.000298201925,0};
  double AMul1Sing[15]   = {6.0131921,-0.451572318,1.06409876,-0.565532819,
    -0.0265785191,0,-0.107782872,-0.166377441,0,0.235821288,0.208159863,
    -0.00831875258,0,-0.000298201925,0.226117001};
  double AMul1Singim[15] = {0,0,1.06409876,0,-0.0265785191,0,0,-0.166377441,0,
    -0.235821288,0,-0.00831875258,0,0.000298201925,0};
  
  double A0Sing[15] = {2.08493611, 0.0595090739,0, 0.026460114, 0, 0,
    0.00906667786, 0, 0, 0, 0.00287090794, 0, 0, 0, 0};
  double A1Sing[15] = {2.08493623, -0.059510451, 0, 0.026460955, 0, 0,
    -0.00906706693, 0, 0, 0, 0.00287106903, 0, 0, 0, 0};
  
  
  double dTA00X[15]   = {8.40366148e-05,1.24371239e-05,-0.00024223864,
    7.32052642e-07,-2.06649326e-05,1.55421909e-08,-1.33035844e-08,
    -8.90141524e-07,3.28549285e-08,-1.99245602e-07,-8.66506018e-09,
    8.76895226e-09,6.77273003e-09,-3.3842657e-08,-1.12365334e-09};
  double dTA00Xim[15] = {0,0,2.28653741e-06,0,2.54344891e-07,-5.70164427e-06,0,
    9.68136519e-09,-4.93159502e-07,1.10213427e-08,0,-1.48696517e-09,
    -1.63296663e-08,3.05594515e-09,-1.47030724e-08};
  double dTA00Y[15]   = {0.00210486595,0.000412647342,2.28653741e-06,
    4.61199423e-05,2.54344891e-07,1.2074341e-05,3.99665825e-06,9.68136519e-09,
    1.64053954e-06,1.53123411e-08,2.88252398e-07,-1.48696517e-09,
    1.63396968e-07,1.89505581e-09,1.2600396e-08};
  double dTA00Yim[15] = {0,0,-0.000278760786,0,-2.68948152e-05,-2.72384346e-07,
    0,-1.62340475e-06,-2.38996886e-08,1.55756718e-08,0,-6.15154603e-08,
    -1.25082386e-09,7.72893677e-09,-3.70533308e-10};
  double dTA00Z[15]   = {0.00762745576,0.00098754963,8.79437465e-06,
    6.75899726e-05,6.39735563e-07,3.11494127e-06,3.23452259e-06,-8.88300647e-09,
    4.6722108e-07,2.31699754e-08,7.01440824e-08,-7.6619317e-09,
    5.01207285e-08,4.29919861e-09,-1.46977779e-08};
  double dTA00Zim[15] = {0,0,0.000291785734,0,3.87091775e-05,2.54344891e-07,0,
    3.54874075e-06,1.26798933e-08,4.68415939e-07,0,2.63294193e-07,
    -2.06841692e-09,7.88415973e-08,1.75044317e-09};
  
  double dTA11X[15]   = {-0.0151853099,0.000252509934,0.00211174413,
    0.000216265654,1.16422283e-05,-0.000293015978,2.59871423e-06,
    -2.97749845e-05,-7.96572085e-07,3.73320959e-05,-3.46401015e-06,
    6.28382238e-09,3.50804631e-06,1.54270798e-08,-4.49835015e-06};
  double dTA11Xim[15] = {0,0,-2.92423424e-05,0,-7.57435756e-06,2.00432626e-05,
    0,1.68562581e-06,-2.72333015e-07,-3.1580328e-06,0,-9.15146217e-08,
    -1.55301981e-07,2.94339673e-08,4.47759661e-07};
  double dTA11Y[15]   = {-0.000649545295,0.000193817791,-2.92423424e-05,
    -2.76416645e-05,-7.57435756e-06,9.08336352e-06,2.0103999e-06,1.68562581e-06,
    8.0280573e-07,-2.27851923e-06,-8.59235284e-08,-9.15146217e-08,
    -2.08106667e-07,-3.86472161e-08,3.7684397e-07};
  double dTA11Yim[15] = {0,0,-0.00116193389,0,1.57819644e-05,0.00019713346,0,
    6.07902481e-06,-5.9292554e-07,-2.75674893e-05,0,5.2794759e-08,
    -1.30545444e-06,1.85252068e-09,3.55102393e-06};
  double dTA11Z[15]   = {-0.00397272143,-0.000901392458,0.000178551487,
    -3.84530678e-05,0.000168987733,-2.06986804e-06,2.69764415e-05,
    2.30411642e-06,-2.23655072e-05,-8.31383539e-08,-7.24390524e-08,
    -3.07335924e-06,-3.31687119e-08,2.55517124e-06,4.79933141e-09};
  double dTA11Zim[15] = {0,0,0.000137049874,0,-2.28940877e-05,-7.57435756e-06,
    0,1.78217115e-06,2.12661856e-06,2.16564579e-07,0,-8.01151965e-08,
    -1.30447713e-07,-1.93195131e-07,-3.25737536e-09};
  
  double dTA22X[15]   = {0.0138345289,-0.000264971459,0.00181443092,
    -0.000171521013,9.31304576e-06,0.000239127813,-2.5849498e-06,
    -2.31391192e-05,7.63301685e-07,2.90615005e-05,2.61981372e-06,
    -1.4916962e-08,-2.62055104e-06,1.90680347e-08,3.33215071e-06};
  double dTA22Xim[15] = {0,0,2.69513189e-05,0,7.82883757e-06,
    1.94376751e-05,0,-1.6952699e-06,-1.71150792e-07,3.1468793e-06,0,
    9.0024276e-08,-1.67087182e-07,-2.63740263e-08,4.4839425e-07};
  double dTA22Y[15]   = {-0.000617451849,0.000184236798,2.69513189e-05,
    -2.60347135e-05,7.82883757e-06,9.11776478e-06,1.80861763e-06,
    -1.6952699e-06,7.94752641e-07,2.26309537e-06,-6.50883031e-08,
    9.0024276e-08,-2.06803856e-07,4.05436927e-08,3.7554146e-07};
  double dTA22Yim[15] = {0,0,-0.00101028362,0,1.51927385e-05,-0.00016218759,0,
    4.55402329e-06,6.16985011e-07,-2.14735717e-05,0,4.94370144e-08,
    9.54333145e-07,4.82955894e-09,-2.62572648e-06};
  double dTA22Z[15]   = {-0.0038492851,-0.000721285509,-0.000187363115,
    -3.47275886e-05,-0.000133129258,-2.93984638e-06,2.11442689e-05,
    -2.29483949e-06,-1.7246713e-05,5.97335319e-08,-4.43405578e-08,
    2.32349965e-06,-4.59035476e-08,-1.89632039e-06,5.03406139e-09};
  double dTA22Zim[15] = {0,0,0.000130275089,0,-2.15659715e-05,7.82883757e-06,0,
    1.6053902e-06,-2.13925215e-06,2.54584389e-07,0,-6.1290081e-08,
    1.28374534e-07,-1.98818893e-07,5.00973358e-09};
  
  double dTA02X[15]   = {-0.00134312154,0.000264445158,0.000116010163,-2.97402876e-05,-9.40166304e-06,6.18828082e-06,2.59867123e-06,3.26725526e-07,-7.75828158e-07,1.28879662e-07,-1.89570995e-07,1.67712122e-08,6.96531629e-08,-2.07034721e-08,-3.0608302e-10};
  double dTA02Xim[15] = {0,0,4.55809629e-05,0,-7.82860105e-06,-2.27082704e-06,0,9.29391901e-07,1.71215149e-07,-1.96522678e-07,0,-9.00289095e-08,-9.13508283e-10,2.63862862e-08,-7.31276472e-09};
  double dTA02Y[15]   = {0.000936010711,-0.000184235667,4.55809629e-05,2.07108202e-05,-7.82860105e-06,5.71633817e-06,-1.80864707e-06,9.29391901e-07,-7.94704561e-07,2.69290902e-07,1.31834659e-07,-9.00289095e-08,8.13075048e-08,-4.05349697e-08,7.61554053e-09};
  double dTA02Yim[15] = {0,0,0.00014957279,0,-1.51711926e-05,5.28159952e-06,0,1.01243297e-06,-6.12044688e-07,4.296989e-08,0,-4.97434464e-08,4.8578427e-08,-4.01876294e-09,-3.60700385e-09};
  double dTA02Z[15]   = {0.00384765899,-0.000501599176,0.000186990964,3.48077424e-05,-2.49543843e-05,2.88476479e-06,-1.71518864e-06,2.30661096e-06,-4.36863404e-07,-6.68643219e-08,4.26242074e-08,-1.73086554e-07,4.7423916e-08,1.129776e-08,-5.89893548e-09};
  double dTA02Zim[15] = {0,0,-0.000130274289,0,1.73783056e-05,-7.82860105e-06,0,-1.60541546e-06,1.18424623e-06,-2.54538487e-07,0,1.20374152e-07,-1.28381054e-07,4.30971946e-08,-5.002315e-09};
  
  double dTA10X[15]   = {-0.00130993449,-0.000254679919,0.000113852059,-2.81182583e-05,9.17038656e-06,6.22042771e-06,-2.39629293e-06,3.19878772e-07,7.83716188e-07,1.13261912e-07,-1.68775884e-07,-1.50020177e-08,7.09628726e-08,1.87642858e-08,-1.60443988e-09};
  double dTA10Xim[15] = {0,0,-4.96111547e-05,0,-8.41478095e-06,2.5486814e-06,0,-9.83275685e-07,1.96040513e-07,2.18428122e-07,0,-9.33933949e-08,2.27461873e-09,2.93631643e-08,7.6774407e-09};
  double dTA10Y[15]   = {-0.00103640083,-0.000201539018,-4.96111547e-05,-2.22575405e-05,-8.41478095e-06,-6.01973614e-06,-1.89755251e-06,-9.83275685e-07,-8.16253293e-07,-2.73962252e-07,-1.33718702e-07,-9.33933949e-08,-8.10510319e-08,-4.00164781e-08,-6.94320742e-09};
  double dTA10Yim[15] = {0,0,0.000137353498,0,1.3152869e-05,4.60385499e-06,0,7.84698074e-07,4.94619605e-07,3.11319409e-09,0,2.90853856e-08,3.4206926e-08,-2.38930338e-09,-5.21757356e-09};
  double dTA10Z[15]   = {-0.00375206577,-0.000479482777,-0.000180085897,-3.19354823e-05,-2.36134099e-05,-1.99124122e-06,-1.43640969e-06,-2.12927475e-06,-2.96207016e-07,1.18023186e-07,-2.101148e-08,-1.54295517e-07,-3.14443593e-08,1.97054319e-08,7.47892317e-09};
  double dTA10Zim[15] = {0,0,-0.000142509607,0,-1.86913911e-05,-8.41478095e-06,0,-1.68608032e-06,-1.25316907e-06,-2.53200807e-07,0,-1.22243709e-07,-1.33218596e-07,-4.22334851e-08,-3.76651522e-09};
  
  double dTA21X[15]   = {0.0139262899,-5.10438415e-07,-0.00224230427,-0.000245273608,8.5965464e-08,0.000299219901,1.33139083e-08,2.91976279e-05,-1.21544176e-08,-3.7433673e-05,3.26576252e-06,-1.79959726e-09,-3.4316127e-06,1.58722952e-09,4.49691933e-06};
  double dTA21Xim[15] = {0,0,7.25351982e-05,0,2.25090129e-10,-1.71675376e-05,0,-7.65908567e-07,-6.12575553e-11,2.9504751e-06,0,-4.41112622e-12,1.68007384e-07,1.16715163e-11,-4.41099212e-07};
  double dTA21Y[15]   = {-0.000318571689,-1.07581356e-09,7.25351982e-05,5.3241063e-06,2.25090129e-10,-1.48346991e-05,2.80247377e-11,-7.65908567e-07,-4.5763443e-11,2.53248806e-06,-6.67490054e-08,-4.41112622e-12,1.25501358e-07,8.30440738e-12,-3.83172417e-07};
  double dTA21Yim[15] = {0,0,0.00100830138,0,-2.09013662e-08,-0.000192124427,0,-7.14553469e-06,4.79360922e-09,2.75463418e-05,0,2.97400437e-10,1.35278166e-06,-7.86897492e-10,-3.55500189e-06};
  double dTA21Z[15]   = {1.57677268e-06,0.00143160944,-3.60934464e-07,-7.7755832e-08,-0.00019330212,5.34334151e-08,-2.48991735e-05,1.14219056e-08,2.26770637e-05,-6.91900298e-09,1.66576792e-09,2.89260044e-06,-1.47555988e-09,-2.53956994e-06,8.39380656e-10};
  double dTA21Zim[15] = {0,0,-7.60715061e-10,0,4.18783338e-06,2.25090129e-10,0,2.40419269e-11,-9.55044034e-07,-4.36911396e-11,0,-5.90864162e-08,-6.20776981e-12,1.55727906e-07,7.06255555e-12};
  
  double dA00[15] = {1.0544637e-06,7.85368589e-05,-3.7739641e-05,1.1710611e-05,
    -4.2581267e-06,-2.54957295e-06,1.1509835e-06,-2.02804679e-07,-3.6652067e-07,
    -4.9837825e-08,8.9749664e-08,2.9501223e-09,-3.6136045e-08,-8.8566072e-09,0};
  double dA01im[15] = {0,0,4.49886716e-05,0,5.90910631e-06,1.96972506e-06,0,
    4.233064e-07,2.48485723e-07,5.45274128e-09,0,1.98818687e-08,1.98691187e-08,
    0,-2.28180158e-09};
  double dA11[15] = {-3.1633912e-06,0.00015605884,5.5533741e-05,-1.438732e-05,
    -9.82898223e-06,-8.25866317e-07,8.00828799e-07,1.0218887e-06,
    1.40543141e-07,-4.82103708e-08,-2.17563271e-08,-8.19726265e-08,
    -1.62881095e-08,8.72542217e-09,3.2623144e-09};
  double dA10im[15] = {0,0,-1.44387654e-05,0,3.28980681e-06,-1.10667758e-06,0,
    -4.39392408e-07,1.06332252e-07,9.822109e-08,0,4.54976027e-08,-3.4394996e-09,
    -1.4438813e-08,-3.54879966e-09};
  
  double dATrip02[15] = {2.5859665e-05,0.000620232366,0.000285436615,
    0.000236529621,8.34048199e-05,-4.91547478e-05,6.02966515e-05,1.00209838e-05,
    -1.8599944e-05,2.4739817e-06,1.2321061e-05,-5.74248967e-07,-4.85496063e-06,
    1.14977575e-06,9.95493175e-08};
  double dATrip11im[15] = {1.26068667e-07,1.77759087e-06,-0.0183902725,
    3.50966837e-06,0.00404630668,0.031857372,-4.92548117e-06,0.00782664265,
    -0.00183388625,-0.0429967544,5.87124061e-06,0.000568727219,-0.018891895,
    -7.10651767e-05,0.051565011};
  double dATrip21[15] = {1.11111691e-06,2.36981299e-05,-2.85249661e-05,
    -2.31735024e-06,2.95909727e-06,2.85184934e-06,-1.8419939e-07,4.31283004e-08,
    1.5601647e-07,3.89577124e-08,4.5967117e-09,-1.62540654e-08,-3.85136191e-09,
    1.48585811e-08,4.41682557e-09};
  double dATrip10im[15] = {-3.62501039e-08,-3.17816832e-07,-0.000808271535,
    -1.97785916e-07,0.00276094976,-0.000692074481,3.07878903e-06,-0.00162795109,
    0.00110762117,0.000646125447,-2.84638226e-06,0.00264295314,0.000193963401,
    -0.000577338053,-0.000156288848};
  
  double dASing00[15] = {0,0,-0.0012172202,0,-0.0009384587,0,0,-0.00045537497,
    0,0,0,-0.000186430392,0,0,0};
  double dASing00im[15] = {0,0,-0.00121722018,0,-0.00093845868,0,0,
    -0.000455374968,0,0,0,-0.000186430392,0,0,0};
  double dASing01[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double dASing01im[15] = {0,0,0.00121722018,0,0.00093845868,0,0,
    0.000455374968,0,0,0,0.000186430392,0,0,0};
  double dASing11[15] = {-0.000411266015,0.0044280306,0,-0.00261408043,0,0,
    0.00113979399,0,0,0,-0.000444103676,0,0,0,0};
  double dASing11im[15] = {0,0,-0.00121722018,0,0.00093845868,0,0,
    -0.000455374968,0,0,0,0.000186430392,0,0,0};

  double L0[15] = {0.0676056845,0.00750146809,0.00185347754,0.000437801593,
    0.000208024957,1.36224131e-05,2.00421651e-05,1.64224844e-05,1.78797695e-06,
    -5.89301269e-07,6.82090982e-07,1.08352568e-06,1.70746019e-07,
    -8.90524479e-08,-2.94499438e-08};
  double L0_im[15] = {0,0,0.0014651663,0,0.000164472954,5.73054024e-05,0,
    1.29876352e-05,7.52793049e-06,1.2692253e-06,0,8.57198339e-07,
    7.19654983e-07,1.91648298e-07,1.49038432e-08};
  double L1[15] = {0.0698331075,-0.00794092114,-0.0017799131,0.000487357134,
    0.000205679316,1.33067114e-05,-2.43860909e-05,-1.69189375e-05,
    -1.80158025e-06,4.36309276e-07,1.00017433e-06,1.17946738e-06,
    1.78774053e-07,-6.80727362e-08,-2.14820881e-08};
  double L1_im[15] = {0,0,-0.00136887958,0,0.00015821673,4.99704248e-05,0,
    -1.30187435e-05,-6.77152753e-06,-1.04925268e-06,0,9.07947072e-07,
    6.72700267e-07,1.63546868e-07,1.23857389e-08};
  
  double L0Sing[15] = {0.0210272529,0.0124511487,0,0.00424544438,0,0,
    0.0013056269,0,0,0,0.000389723882,0,0,0,0};
  double L0SingIm[15] = {};
  double L1Sing[15] = {0.0210272529,-0.0124511487,0,0.00424544438,0,0,
    -0.0013056269,0,0,0,0.000389723882,0,0,0,0};
  double L1SingIm[15] = {};
  
  double dLdx0[15] = {1.6750432e-06,-6.48184643e-05,2.34002566e-05,
    9.24707607e-06,1.70674585e-05,4.63747276e-06,2.14310222e-06,4.22726785e-06,
    1.31460586e-06,1.10177468e-07,-4.90918167e-07,-2.74985149e-06,
    -8.99415319e-07,1.44385542e-07,1.00536605e-07};
  double dAdy0im[15] = {-1.61770896e-06,-3.41962773e-06,3.5234244e-05,
    -1.05669019e-05,7.2297967e-06,4.09008335e-06,-3.47883114e-06,7.54460728e-07,
    1.96413284e-06,4.91172588e-07,-1.26891901e-08,-2.75940448e-06,
    -1.43810291e-06,-9.21941465e-08,6.53481772e-08};
  double dLdz0[15] = {-0.000150691983,-1.33519347e-05,-2.06166021e-05,
    5.03159396e-05,2.06163519e-05,1.07343682e-06,2.0450875e-05,1.79780738e-05,
    2.12705147e-06,-5.04752566e-07,-4.9965286e-06,-4.07149703e-06,
    -6.09290158e-07,1.72960474e-07,5.43921381e-08};
  double dLdx1[15] = {-0.000987674279,0.000187867261,-0.000357770967,
    -0.000293354883,5.61544215e-05,4.6981816e-05,3.09724418e-05,-3.49290618e-05,
    -2.36043496e-05,-2.24351551e-06,-2.44654805e-05,-6.57150526e-06,
    5.77399697e-06,1.98995506e-06,-1.03763066e-08};
  double dAdy1im[15] = {2.61723492e-07,-3.77632927e-08,-0.00037629043,
    9.4750488e-08,9.7902191e-05,3.25167845e-05,-1.47164999e-08,-3.90081e-05,
    -2.15522995e-05,-1.42788374e-06,1.30984974e-08,1.38167558e-07,
    9.21538861e-07,-7.45070972e-07,-4.6128438e-07};
  double dAdz1im[15] = {-1.18032719e-07,-7.74616576e-08,0.000240731667,
    7.46244492e-08,-0.000116771598,-5.26094072e-05,-4.6568393e-08,
    4.78806389e-05,3.23462883e-05,6.20808799e-06,1.75724926e-08,-7.2420955e-06,
    -8.11121113e-06,-2.60581756e-06,-2.39864327e-07};
  
  double dLdx0Sing[15] = {0,0,-0.00025467424,0,-0.000150567895,0,0,
    -6.5572367e-05,0,0,0,-2.53063413e-05,0,0,0};
  double dAdy0imSing[15] = {0,0,-0.00025467424,0,-0.000150567895,0,0,
    -6.5572367e-05,0,0,0,-2.53063413e-05,0,0,0};
  double dLdz0Sing[15] = {0.00128088088,0.000940448055,0,0.000393219492,0,0,
    0.000151918661,0,0,0,5.43212252e-05,0,0,0,0};
  double dLdx1Sing[15] = {0,0,-0.000254674238,0,0.000150567895,0,0,
    -6.5572367e-05,0,0,0,2.53063413e-05,0,0,0};
  double dAdy1imSing[15] = {0,0,-0.00025467424,0,0.00015056789,0,0,
    -6.5572367e-05,0,0,0,2.53063413e-05,0,0,0};
  double dAdz1imSing[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

} ; // end ASolverUTest


TEST_F(ASolverUTest, checkGamma)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);

  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 1).real(),  1.463995711, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 0, 5).real(),  1.760111936, preclim);
  
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 2).real(),  1.621243794, preclim);
  EXPECT_NEAR( ASolvTest.get_gamma_ni( 1, 7).real(),  1.799701878, preclim);
}

TEST_F(ASolverUTest, checkDelta)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 1).real(),   0.87554313, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 0, 5).real(),   0.06832297, preclim);
  
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 2).real(),   11.4370663, preclim);
  EXPECT_NEAR( ASolvTest.get_delta_ni( 1, 7).real()/181.9847, 1.0, preclim);
}


TEST_F(ASolverUTest, checkE)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  
  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 0, 0).real(), 5.0, preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 5, 0).real()/-4.7683716e-06, 1.0, preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 5, 0).imag(),                0.0, preclim);

  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 6, -5).real(),        0.0, preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 0, 6, -5).imag(),        0.0, preclim);

  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 0, 0).real(),-0.4, preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 3, -2).real(), 0.00014227949, preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 3, -2).imag(), 0.00134791098, preclim);

  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 6, -5).real()/-6.7122985e-06, 1., preclim);
  EXPECT_NEAR(ASolvTest.get_E_ni( 1, 6, -5).imag()/-3.8653666e-06, 1., preclim);
}

TEST_F(ASolverUTest, checkSH)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 0, 0).real(), 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 5, 0).real(),-1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 5, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 6,-5).real(), 0.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 0, 0, 6,-5).imag(), 0.0, preclim);
  
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 0, 0).real(), 1.0, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 0, 0).imag(), 0.0, preclim);

  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 3,-2).real(),-0.0522110883, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 3,-2).imag(), 0.4946303982, preclim);
  
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 6, 5).real(), 0.3615486465, preclim);
  EXPECT_NEAR( ASolvTest.get_SH_ij( 1, 0, 6, 5).imag(),-0.2082023636, preclim);
}

TEST_F(ASolverUTest, checkAMulti)
{
  mol_.clear( );
  Pt pos[3] = { Pt( 0.0, 0.0, -5.0 ),
    Pt( 10.0, 7.8, 25.0 ),Pt(-10.0, 7.8, 25.0) };
  for (int molInd = 0; molInd < 3; molInd ++ )
  {
    int M = 3;
    vector<double> charges(M); vector<Pt> posCharges(M);
    charges[0]    = 2.0; posCharges[0] = pos[molInd];
    charges[1]    = 2.0; posCharges[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]    = 2.0; posCharges[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
    mol_.push_back( molNew );
  }
  
  const int vals           = nvals;
  int nmol                 = 3;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20);

  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR(ASolvTest.get_A_ni(0, n, m).real()/ATrip0[ct], 1.0, preclim);
      EXPECT_NEAR(ASolvTest.get_A_ni(1, n, m).real()/ATrip1[ct], 1.0, preclim);
      EXPECT_NEAR(ASolvTest.get_A_ni(2, n, m).real()/ATrip2[ct], 1.0, preclim);
      
      if (ATrip0im[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(0, n, m).imag()/ATrip0im[ct],
                    1.0, preclim);
      if (ATrip1im[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(1, n, m).imag()/ATrip1im[ct],
                    1.0, preclim);
      if (ATrip2im[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(2, n, m).imag()/ATrip2im[ct],
                    1.0, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkA)
{
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20);

  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR(ASolvTest.get_A_ni( 0, n, m).real()/A0[ct],      1, preclim);
      if (A0_im[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni( 0, n, m).imag()/A0_im[ct], 1, preclim);
      EXPECT_NEAR(ASolvTest.get_A_ni( 1, n, m).real()/A1[ct],      1, preclim);
      if (A1_im[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni( 1, n, m).imag()/A1_im[ct], 1, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkASingMult)
{
  mol_sing_.clear( );
  Pt pos[3] = { Pt( 0.0, 0.0, -5.0 ), Pt( 0.0, 0.0, 0.0 )};
  for (int molInd = 0; molInd < 2; molInd ++ )
  {
    int M = 3; vector<double> charges(M); vector<Pt> posCharges(M);
    charges[0]    = 2.0; posCharges[0] = pos[molInd];
    charges[1]    = 2.0; posCharges[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
    charges[2]    = 2.0; posCharges[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
    
    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
    mol_sing_.push_back( molNew );
  }
  
  const int vals           = 5;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-20);

  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      if (AMul0Sing[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(0,n,m).real()/AMul0Sing[ct], 1, preclim);
      if (AMul0Singim[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(0,n,m).imag()/AMul0Singim[ct],
                    1, preclim);
      if (AMul1Sing[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(1,n,m).real()/AMul1Sing[ct], 1, preclim);
      if (AMul1Singim[ct] != 0)
        EXPECT_NEAR(ASolvTest.get_A_ni(1,n,m).imag()/AMul1Singim[ct],
                    1, preclim);
      ct++;
    }
  }
}

TEST_F(ASolverUTest, checkASing)
{
  const int vals           = 5;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-30);
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).real(), A0Sing[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 0, n, m).imag(),        0.0, preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).real(), A1Sing[ct], preclim);
      EXPECT_NEAR( ASolvTest.get_A_ni( 1, n, m).imag(),        0.0, preclim);
      ct++;
    }
  }
}

//TEST_F(ASolverUTest, checkgradT_A)
//{
//  mol_.clear( );
//  Pt pos[3] = { Pt( 0.0, 0.0, -5.0 ),
//    Pt( 10.0, 7.8, 25.0 ),Pt(-10.0, 7.8, 25.0) };
//  for (int molInd = 0; molInd < 3; molInd ++ )
//  {
//    int M = 3;
//    vector<double> charges(M); vector<Pt> posCharges(M);
//    charges[0]    = 2.0; posCharges[0] = pos[molInd];
//    charges[1]    = 2.0; posCharges[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
//    charges[2]    = 2.0; posCharges[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
//    
//    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
//    mol_.push_back( molNew );
//  }
//  
//  const int vals           = 10;
//  int nmol                 = 3;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-26);
//  
//  VecOfMats<cmplx>::type dT_A00 = ASolvTest.get_gradT_Aij( 0, 0);
//  VecOfMats<cmplx>::type dT_A11 = ASolvTest.get_gradT_Aij( 1, 1);
//  VecOfMats<cmplx>::type dT_A22 = ASolvTest.get_gradT_Aij( 2, 2);
//  
//  VecOfMats<cmplx>::type dT_A02 = ASolvTest.get_gradT_Aij( 0, 2);
//  VecOfMats<cmplx>::type dT_A10 = ASolvTest.get_gradT_Aij( 1, 0);
//  VecOfMats<cmplx>::type dT_A21 = ASolvTest.get_gradT_Aij( 2, 1);
//  
//  int ct = 0;
//  for ( int n = 0; n < 5; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
//      if (dTA00X[ct] != 0) EXPECT_NEAR(dT_A00[0](n, m+nvals).real()/dTA00X[ct],
//                                       1.0, preclim);
//      if (dTA00Y[ct] != 0) EXPECT_NEAR(dT_A00[1](n, m+nvals).real()/dTA00Y[ct],
//                                       1.0, preclim);
//      if (dTA00Z[ct] != 0) EXPECT_NEAR(dT_A00[2](n, m+nvals).real()/dTA00Z[ct],
//                                       1.0, preclim);
//      if (dTA00Xim[ct] != 0) EXPECT_NEAR(dT_A00[0](n, m+nvals).imag()/
//                                         dTA00Xim[ct], 1.0, preclim);
//      if (dTA00Yim[ct] != 0) EXPECT_NEAR(dT_A00[1](n, m+nvals).imag()/
//                                         dTA00Yim[ct], 1.0, preclim);
//      if (dTA00Zim[ct] != 0) EXPECT_NEAR(dT_A00[2](n, m+nvals).imag()/
//                                         dTA00Zim[ct], 1.0, preclim);
//      
//      if (dTA11X[ct] != 0) EXPECT_NEAR(dT_A11[0](n, m+nvals).real()/dTA11X[ct],
//                                       1.0, preclim);
//      if (dTA11Y[ct] != 0) EXPECT_NEAR(dT_A11[1](n, m+nvals).real()/dTA11Y[ct],
//                                       1.0, preclim);
//      if (dTA11Z[ct] != 0) EXPECT_NEAR(dT_A11[2](n, m+nvals).real()/dTA11Z[ct],
//                                       1.0, preclim);
//      if (dTA11Xim[ct] != 0) EXPECT_NEAR(dT_A11[0](n, m+nvals).imag()/
//                                         dTA11Xim[ct], 1.0, preclim);
//      if (dTA11Yim[ct] != 0) EXPECT_NEAR(dT_A11[1](n, m+nvals).imag()/
//                                         dTA11Yim[ct], 1.0, preclim);
//      if (dTA11Zim[ct] != 0) EXPECT_NEAR(dT_A11[2](n, m+nvals).imag()/
//                                         dTA11Zim[ct], 1.0, preclim);
//      
//      if (dTA22X[ct] != 0) EXPECT_NEAR(dT_A22[0](n, m+nvals).real()/dTA22X[ct],
//                                       1.0, preclim);
//      if (dTA22Y[ct] != 0) EXPECT_NEAR(dT_A22[1](n, m+nvals).real()/dTA22Y[ct],
//                                       1.0, preclim);
//      if (dTA22Z[ct] != 0) EXPECT_NEAR(dT_A22[2](n, m+nvals).real()/dTA22Z[ct],
//                                       1.0, preclim);
//      if (dTA22Xim[ct] != 0) EXPECT_NEAR(dT_A22[0](n, m+nvals).imag()/
//                                         dTA22Xim[ct], 1.0, preclim);
//      if (dTA22Yim[ct] != 0) EXPECT_NEAR(dT_A22[1](n, m+nvals).imag()/
//                                         dTA22Yim[ct], 1.0, preclim);
//      if (dTA22Zim[ct] != 0) EXPECT_NEAR(dT_A22[2](n, m+nvals).imag()/
//                                         dTA22Zim[ct], 1.0, preclim);
//      
//      if (dTA02X[ct] != 0) EXPECT_NEAR(dT_A02[0](n, m+nvals).real()/dTA02X[ct],
//                                       1.0, preclim);
//      if (dTA02Y[ct] != 0) EXPECT_NEAR(dT_A02[1](n, m+nvals).real()/dTA02Y[ct],
//                                       1.0, preclim);
//      if (dTA02Z[ct] != 0) EXPECT_NEAR(dT_A02[2](n, m+nvals).real()/dTA02Z[ct],
//                                       1.0, preclim);
//      if (dTA02Xim[ct] != 0) EXPECT_NEAR(dT_A02[0](n, m+nvals).imag()/
//                                         dTA02Xim[ct], 1.0, preclim);
//      if (dTA02Yim[ct] != 0) EXPECT_NEAR(dT_A02[1](n, m+nvals).imag()/
//                                         dTA02Yim[ct], 1.0, preclim);
//      if (dTA02Zim[ct] != 0) EXPECT_NEAR(dT_A02[2](n, m+nvals).imag()/
//                                         dTA02Zim[ct], 1.0, preclim);
//      
//      if (dTA10X[ct] != 0) EXPECT_NEAR(dT_A10[0](n, m+nvals).real()/dTA10X[ct],
//                                       1.0, preclim);
//      if (dTA10Y[ct] != 0) EXPECT_NEAR(dT_A10[1](n, m+nvals).real()/dTA10Y[ct],
//                                       1.0, preclim);
//      if (dTA10Z[ct] != 0) EXPECT_NEAR(dT_A10[2](n, m+nvals).real()/dTA10Z[ct],
//                                       1.0, preclim);
//      if (dTA10Xim[ct] != 0) EXPECT_NEAR(dT_A10[0](n, m+nvals).imag()/
//                                         dTA10Xim[ct], 1.0, preclim);
//      if (dTA10Yim[ct] != 0) EXPECT_NEAR(dT_A10[1](n, m+nvals).imag()/
//                                         dTA10Yim[ct], 1.0, preclim);
//      if (dTA10Zim[ct] != 0) EXPECT_NEAR(dT_A10[2](n, m+nvals).imag()/
//                                         dTA10Zim[ct], 1.0, preclim);
//      
//      if (dTA21X[ct] != 0) EXPECT_NEAR(dT_A21[0](n, m+nvals).real()/dTA21X[ct],
//                                       1.0, preclim);
//      if (dTA21Y[ct] != 0) EXPECT_NEAR(dT_A21[1](n, m+nvals).real()/dTA21Y[ct],
//                                       1.0, preclim);
//      if (dTA21Z[ct] != 0) EXPECT_NEAR(dT_A21[2](n, m+nvals).real()/dTA21Z[ct],
//                                       1.0, preclim);
//      if (dTA21Xim[ct] != 0) EXPECT_NEAR(dT_A21[0](n, m+nvals).imag()/
//                                         dTA21Xim[ct], 1.0, preclim);
//      if (dTA21Yim[ct] != 0) EXPECT_NEAR(dT_A21[1](n, m+nvals).imag()/
//                                         dTA21Yim[ct], 1.0, preclim);
//      if (dTA21Zim[ct] != 0) EXPECT_NEAR(dT_A21[2](n, m+nvals).imag()/
//                                         dTA21Zim[ct], 1.0, preclim);
//      ct++;
//    }
//  }
//}


//TEST_F(ASolverUTest, checkgradA)
//{
//  mol_.clear( );
//  Pt pos[3] = { Pt( 0.0, 0.0, -5.0 ),
//    Pt( 10.0, 7.8, 25.0 ),Pt(-10.0, 7.8, 25.0) };
//  for (int molInd = 0; molInd < 3; molInd ++ )
//  {
//    int M = 3;
//    vector<double> charges(M); vector<Pt> posCharges(M);
//    charges[0]    = 2.0; posCharges[0] = pos[molInd];
//    charges[1]    = 2.0; posCharges[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
//    charges[2]    = 2.0; posCharges[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
//    
//    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
//    mol_.push_back( molNew );
//  }
//  
//  const int vals           = 5;
//  int nmol                 = 3;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-15); ASolvTest.solve_gradA(1E-36);
//  
////  for ( int n = 0; n < nmol; n++ )
////  {
////    for ( int m = 0; m < nmol; m++ )
////    {
////      ASolvTest.print_dAidx(m, m, 5);
////      ASolvTest.print_dAidy(m, m, 5);
////      ASolvTest.print_dAidz(m, m, 5);
////    }
////  }
//  
//  int ct = 0;
//  for ( int n = 0; n < 3; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
////      if (dATrip02[ct] != 0)
////        EXPECT_NEAR( ASolvTest.get_dAdx_ni( 0, 2, n, m).real()/dATrip02[ct],
////                     1.0, preclim);
////      if (dATrip11im[ct] != 0)
////        EXPECT_NEAR( ASolvTest.get_dAdy_ni( 1, 1, n, m).imag()/dATrip11im[ct],
////                     1.0, preclim);
////      if (dATrip21[ct] != 0)
////        EXPECT_NEAR( ASolvTest.get_dAdz_ni( 2, 1, n, m).real()/dATrip21[ct],
////                     1.0, preclim);
////      if (dATrip10im[ct] != 0)
////        EXPECT_NEAR( ASolvTest.get_dAdx_ni( 1, 0, n, m).imag()/dATrip10im[ct],
////                     1.0, preclim);
//      ct++;
//    }
//  }
//}

TEST_F(ASolverUTest, checkgradASing)
{
//  mol_sing_.clear();
//  Pt cgPosSi[2] = { Pt( 0.0, 0.0, -5.0 ), Pt( 0.0, 0.0, 0.0 ) };
//  
//  for (int molInd = 0; molInd < 2; molInd ++ )
//  {
//    int M = 1; vector<double> charges(1); vector<Pt> posCharges(1);
//    
//    charges[0]    = 2.0;
//    posCharges[0] = cgPosSi[molInd];
//    
//    Molecule molSing( M, 2.0, charges, posCharges);
//    mol_sing_.push_back( molSing );
//  }
  
  const int vals           = nvals;
  int nmol                 = 2;
  BesselConstants bConsta  = BesselConstants( 2*vals );
  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
  System sys               = System( const_, mol_sing_ );
  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
                                      sys.get_lambda(), nvals);
  
  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
  ASolvTest.solve_A(1E-30); ASolvTest.solve_gradA(1E-16);
  
//  for ( int m = 0; m < nmol; m++ )
//  {
//    ASolvTest.print_dAidx(m, m, 5);
//    ASolvTest.print_dAidy(m, m, 5);
//    ASolvTest.print_dAidz(m, m, 5);
//  }
  
  int ct = 0;
  for ( int n = 0; n < 5; n++ )
  {
    for ( int m = 0; m <= n; m++ )
    {
//      EXPECT_NEAR( ASolvTest.get_dAdx_ni( 0, 0, n, m).real(),
//                  dASing00[ct], preclim);
//      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 0, n, m).imag(),
//                  dASing00im[ct], preclim);
//      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 1, n, m).real(),
//                  dASing01[ct], preclim);
//      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 0, 1, n, m).imag(),
//                  dASing01im[ct], preclim);
//      EXPECT_NEAR( ASolvTest.get_dAdz_ni( 1, 1, n, m).real(),
//                  dASing11[ct], preclim);
//      EXPECT_NEAR( ASolvTest.get_dAdy_ni( 1, 1, n, m).imag(),
//                  dASing11im[ct], preclim);
      ct++;
    }
  }
}

//TEST_F(ASolverUTest, checkL)
//{
//  mol_.clear( );
//  Pt pos[3] = { Pt( 0.0, 0.0, -5.0 ),
//    Pt( 10.0, 7.8, 25.0 ),Pt(-10.0, 7.8, 25.0) };
//  for (int molInd = 0; molInd < 2; molInd ++ )
//  {
//    int M = 3;
//    vector<double> charges(M); vector<Pt> posCharges(M);
//    charges[0]    = 2.0; posCharges[0] = pos[molInd];
//    charges[1]    = 2.0; posCharges[1] = pos[molInd] + Pt(1.0, 0.0, 0.0);
//    charges[2]    = 2.0; posCharges[2] = pos[molInd] + Pt(0.0, 1.0, 0.0);
//    
//    Molecule molNew( M, 2.0, charges, posCharges, pos[molInd]);
//    mol_.push_back( molNew );
//  }
//  const int vals           = nvals;
//  int nmol                 = 2;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-20);
//  VecOfMats<cmplx>::type myL = ASolvTest.calc_L();
//  
//  int ct = 0;
//  for ( int n = 0; n < 5; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
//      EXPECT_NEAR( myL[0](n,m+nvals).real()/L0[ct],    1.0, preclim);
//      if (L0_im[ct] != 0)
//        EXPECT_NEAR( myL[0](n,m+nvals).imag()/L0_im[ct], 1.0, preclim);
//      EXPECT_NEAR( myL[1](n,m+nvals).real()/L1[ct],    1.0, preclim);
//      if (L1_im[ct] != 0)
//        EXPECT_NEAR( myL[1](n,m+nvals).imag()/L1_im[ct], 1.0, preclim);
//      ct++;
//    }
//  }
//}
//
//TEST_F(ASolverUTest, checkLSing)
//{
//  const int vals           = nvals;
//  int nmol                 = 2;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_sing_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-20);
//  
//  VecOfMats<cmplx>::type myL = ASolvTest.calc_L();
//  
//  int ct = 0;
//  for ( int n = 0; n < 5; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
//      EXPECT_NEAR( myL[0](n,m+nvals).real(), L0Sing[ct], preclim);
//      EXPECT_NEAR( myL[0](n,m+nvals).imag(),          0, preclim);
//      EXPECT_NEAR( myL[1](n,m+nvals).real(), L1Sing[ct], preclim);
//      EXPECT_NEAR( myL[1](n,m+nvals).imag(),          0, preclim);
//      ct++;
//    }
//  }
//}

//TEST_F(ASolverUTest, checkdL)
//{
//  const int vals           = nvals;
//  int nmol                 = 2;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
//  MyVector<VecOfMats<cmplx>::type > mydL = ASolvTest.calc_gradL();
//  
//  int ct = 0;
//  for ( int n = 0; n < 5; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
//      if (dLdx0[ct] != 0)
//        EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0[ct],   preclim);
//      if (dAdy0im[ct] != 0)
//        EXPECT_NEAR( mydL[0][1](n,m+nvals).imag(), dAdy0im[ct], preclim);
//      if (dLdx0[ct] != 0)
//        EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0[ct],   preclim);
//      if (dLdx1[ct] != 0)
//        EXPECT_NEAR( mydL[1][0](n,m+nvals).real(), dLdx1[ct],   preclim);
//      if (dAdy1im[ct] != 0)
//        EXPECT_NEAR( mydL[1][1](n,m+nvals).imag(), dAdy1im[ct], preclim);
//      if (dAdz1im[ct] != 0)
//        EXPECT_NEAR( mydL[1][2](n,m+nvals).imag(), dAdz1im[ct], preclim);
//      ct++;
//    }
//  }  
//}

//TEST_F(ASolverUTest, checkdLSing)
//{
//  const int vals           = nvals;
//  int nmol                 = 2;
//  BesselConstants bConsta  = BesselConstants( 2*vals );
//  BesselCalc bCalcu        = BesselCalc( 2*vals, bConsta );
//  SHCalcConstants SHConsta = SHCalcConstants( 2*vals );
//  SHCalc SHCalcu           = SHCalc( 2*vals, SHConsta );
//  System sys               = System( const_, mol_sing_ );
//  ReExpCoeffsConstants re_exp_consts (sys.get_consts().get_kappa(),
//                                      sys.get_lambda(), nvals);
//  
//  ASolver ASolvTest        = ASolver( nmol, vals, bCalcu, SHCalcu, sys);
//  ASolvTest.solve_A(1E-20); ASolvTest.solve_gradA(1E-20);
//  
//  MyVector<VecOfMats<cmplx>::type > mydL = ASolvTest.calc_gradL();
//  
//  int ct = 0;
//  for ( int n = 0; n < 5; n++ )
//  {
//    for ( int m = 0; m <= n; m++ )
//    {
//      EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0Sing[ct],   preclim);
//      EXPECT_NEAR( mydL[0][1](n,m+nvals).imag(), dAdy0imSing[ct], preclim);
//      EXPECT_NEAR( mydL[0][0](n,m+nvals).real(), dLdx0Sing[ct],   preclim);
//      EXPECT_NEAR( mydL[1][0](n,m+nvals).real(), dLdx1Sing[ct],   preclim);
//      EXPECT_NEAR( mydL[1][1](n,m+nvals).imag(), dAdy1imSing[ct], preclim);
//      EXPECT_NEAR( mydL[1][2](n,m+nvals).imag(), dAdz1imSing[ct], preclim);
//      ct++;
//    }
//  }
//}
//

#endif
