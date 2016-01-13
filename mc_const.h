#ifndef _MC_CONST_H
#define _MC_CONST_H 1

//--- PHYSICAL CONSTANTS and CONVERSION FACTORS -----------

const double BOHRRADIUS = 0.5291772108;       // angstrom
const double HARTREE2JL = 4.359748e-18;    // hartree to joule  conversion factor
const double HARTREE2KL = 3.157732e+05;    // hartree to Kelvin conversion factor
const double CMRECIP2KL = 1.4387672;       // cm^-1 to Kelvin conversion factor
const double MHZ2RCM    = 3.335640952e-5;  // MHz to cm^-1 conversion factor

const double HBAR  = 1.05457266; //  (10^-34 Js)     Planck constant
const double AMU   = 1.6605402;  //  (10^-27 kg)     atomic mass unit 
const double K_B   = 1.380658;   //  (10^-23 JK^-1)  Boltzmann constant 
const double WNO2K = 0.6950356; // conversion from CM-1 to K

//--------------- MASSES-----------------------------------

// NIST  [http://www.physics.nist.gov/]
//  
// He3    3.016 029 309 7(9) 
// He4    4.002 603 2497(10)
 
// H      1.007 825 032 1(4) 	  
// D      2.014 101 778 0(4) 	  
// T      3.016 049 2675(11

// C12    12.000 000 0(0) 
// C13    13.003 354 8378(10) 
// C14    14.003 241 988(4)
 
// N14    14.003 074 005 2(9) 
// N15    15.000 108 898 4(9)

// O16    15.994 914 6221(15) 
// O17    16.999 131 50(22) 
// O18    17.999 160 4(9)

// S32    31.972 070 69(12) 
// S33    32.971 458 50(12) 
// S34    33.967 866 83(11) 
// S36    35.967 080 88(25)  

const double MASS_H1   = 1.0078;   // amu
const double MASS_H2   = 2.015650642;   // amu
const double MASS_HE4  = 4.0026032497;   // amu
const double MASS_C12  = 12.0;     // amu
const double MASS_N14  = 14.003;   // amu
const double MASS_O16  = 15.994915;   // amu
const double MASS_S32  = 31.972;   // amu

//----------- ROTATIONAL CONSTANTS --------------------------

const double B_N2O   = 0.602861;  // [K]  N2O  = 0.4190098 cm^-1
const double B_OCS   = 0.292;     // [K]  OCS  
//const double B_CO2   = 0.561122;  // [K]  CO2  = 0.3900 cm^-1 xie  excited(v3=1) 
//const double B_CO2   = 0.557009;  // [K]  CO2  = 0.387141 cm^-1 exp. excited(v3=1) 
const double B_CO2   = 0.561437;  // [K]  CO2  = 0.390219 cm^-1   exp. ground(v3=0) 
const double B_CO   = 2.766086924;  // [K]  CO  = 1.922528955 cm^-1   exp. ground(v=0) 
const double B_HCN   = 2.12682;   // [K]  1.478 221 834 cm^-1 JCP 114 851 (2001)

const double B_HCCCN = 0.218317;  // [K]  4549 MHz [JCP 119 8379 (2003)] approx 0.1517383 cm^-1

// known atom/molecule types

const char HE4[]   =  "He4";
const char H2[]   =  "H2";
const char OCS[]   =  "OCS";
const char N2O[]   =  "N2O";
const char CO2[]   =  "CO2";
const char CO[]   =  "CO";
const char HCN[]   =  "HCN";

const char HCCCN[] =  "HCCCN";
const char H2O[] = "H2O";
const char SO2[] = "SO2";
const char HCOOCH3[] = "HCOOCH3";

#endif  //MC_const.h
