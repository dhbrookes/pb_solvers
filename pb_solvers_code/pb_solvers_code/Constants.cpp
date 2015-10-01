//
//  Constants.cpp
//  pbam
//
//  Created by David Brookes on 9/22/15.
//  Copyright © 2015 David Brookes. All rights reserved.
//

#include "Constants.h"

const double Constants::PERMITTIVITY_VAC = 8.854187817e-12;  // [F/m]=[C2/J/m] , F = Farad
const double Constants::KB = 1.380658e-23;  //!<  [ m^2 kg/ s^2 / K ] = [ J/K ]
const double Constants::LITRE = 1e-3;  // [ m^3/L]
const double Constants::PI = 3.141592654;
const double Constants::COLOUMB_CONSTANT = 8.988e9;  //!< [ N*m^2/C^2 ]
const double Constants::ELECTRON_CHARGE = 1.60217733e-19;  //!<  [ coulombs ]
const double Constants::E2 = 1.60217733e-19 * 1.60217733e-19;
const double Constants::AVOGADRO_NUM = 6.02209e23;  //[C2]
const double Constants::KCAL = 4184.0;  //!<  [ 1 kCal = 4184 Joules ]
const double Constants::ANGSTROM = 1e-10;  //!<  [ 1A = 1e-10 Meters ]
const double Constants::PICO_SEC = 1e-12;  //!<  [ 1 ps = 1e-12 s ]


/*
 Constructor sets default values of independent constants
 */
Constants::Constants()
:bDist_(100.0), qDist_(500.0), fDist_(100.0), dielectricWater_(78.0),
dielectricProt_(4.0), saltConcentration_(0.0100), temp_(353.0), tol_(2.5),
patchAngle_(6.0), rotateAngle_(20.0)
{
    update_all();
}


void Constants::update_kbt()
{
    KbT_ = KB * AVOGADRO_NUM * temp_ / KCAL;
    iKbT_ = 1 / KbT_;
}

void Constants::update_kappa()
{
    double kap_num = sqrt(2 * saltConcentration_ * AVOGADRO_NUM * E2);
    double kap_den = sqrt(LITRE * dielectricWater_ * PERMITTIVITY_VAC * KB * temp_);
    
    kappa_ = ANGSTROM * kap_num * kap_den;
}

void Constants::update_patch_size()
{
    patchSize_ = cos((patchAngle_ * PI) / 180.0);
}

void Constants::update_rotate_size()
{
    patchSize_ = cos((rotateAngle_ * PI) / 180.0);
}

void Constants::update_all()
{
    update_kappa();
    update_rotate_size();
    update_kbt();
    update_rotate_size();
}

const double Constants::convert_int_to_kcal_mol(double val)
{
    double coul_num = E2 * AVOGADRO_NUM;
    double coul_den = PERMITTIVITY_VAC * 4.0 * PI * ANGSTROM * KCAL;
    return val * (coul_num / coul_den );
}

const double Constants::convert_int_to_jmol(double val)
{
    double coul_num = E2 * AVOGADRO_NUM;
    double intj_den = PERMITTIVITY_VAC * 4.0 * PI * ANGSTROM; //IU units density
    return val * (coul_num / intj_den);
}
                                                
                                                
                                                
                                                
