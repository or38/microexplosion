/******************************************************************
This is code for simulation of water in n-dodecane droplet heating
to predict times to puffing

Copyright (c) 2018 Oyuna Rybdylova - All Rights Reserved
You may use, distribute and modify this code under the terms of the MIT license
Please see LICENSE file in the repository https://github.com/or38/microexplosion

If you use the software, please cite
https://github.com/or38/microexplosion
******************************************************************/

#include <math.h>
#include <iostream> 
#include <sstream>
#include <vector>
#include <cstring>

#define PI 3.1415926535897932384626433832795
#define R_air 287.01625988193461525183829875375
#define R_UNIVERSAL 8.31451
#define ACCURACY 1.e-9



using namespace std;
//general subroutines
int Ini();
int Heat_Mass();
int find_Lambda();
int New_thickness();
int Report_T(double ttime);
int Report_main(FILE * ff);
double Update_Ts();


//physical properties subroutines
double Viscosity_gas(double T);
double T_Conductivity_gas(double T);
double Heat_Capacity_gas(double T);

//double Diesel_liquid_density(double T);
//double Diesel_liquid_thermcond(double T);
//double Diesel_liquid_Heat_Capacity(double T);
double Saturation_Pressure(double T);
double kappa0, rad0, Fo_w;
double Evaporation_Heat(double T);
double Diffusion_coef_Vap(double T, double p);
double Evaporation_Heat_ic8(double T);
double Evaporation_Heat_3mp(double T);
double Evaporation_Heat_c7(double T);
double Evaporation_Heat_c16(double T);
double liquid_density(double T);
//double liquid_thermcond(double T);
//double liquid_Heat_Capacity(double T);
double Liquid_diffusivity(double T);
double Saturation_Pressure_ic8(double T);
double Saturation_Pressure_3mp(double T);
double Saturation_Pressure_c7(double T);
double Saturation_Pressure_c16(double T);
double Saturation_Pressure_h2o(double T);
double Saturation_Pressure_diesel(double T);
double Evaporation_Heat_h2o(double T);
double Evaporation_Heat_ndodecane(double T);
double vapour_heatcap(double T);

double density_water(double T);
double heat_capacity_water(double T);
double thermal_conductivity_water(double T);

double density_diesel(double T);
double heat_capacity_diesel(double T);
double thermal_conductivity_diesel(double T);

//double vapour_heatcap(double T);

int Update_vol_fracs();
int Update_mass_fracs();
int print_lambda();



//general variables
char run_name[20];
double time_step, Fo_ts;
double time = 0.0e-15;
double time_end;
int N_Lambda;
int N_Lambda_y;
int N_INT;
int N_thresh;
double Delta_Rw;
double Delta_Rf;
int iprint;
double check_A;
double check_khi;
double check_Tc, check_Tr, check_chi;


//physical variables
double T_gas;
double T_ini;
double P_gas;
double Thickness_ini;
double Thickness;
double delta_dot;
double h0;
double T_eff;
double BM, BT, tot_vap_rate;
double Lewis;
double Sh_Star, Nu_star, Nu, Sh;
double water_vol_ini, vol_w, vol_f;
double water_Rad;
double Y_w, Y_f, x_vs, Ys;
double ratio;
double boundary1 = 1.0;
double boundary2 = 1.0;

vector <double> Temperature_water(1);
vector <double> Temperature_fuel(1);
vector <double> Temperature_water1(1);
vector <double> Temperature_fuel1(1);
double T_w_av, T_f_av;
double T_ref;
double boundary_cond1 = 0.0;
double boundary_cond2 = 0.0;
double dens_old, dens_new;

double evaporation_mass = 0.e-15;
double delta_R = 0.0e-15;
double droplet_mass = 0.e-15;

vector <double> lambda(1);
vector <double> series(1);


//n-dodecane
double c12_mw = 170.33;
double c8_mw = 114.23;
double T_cr_noctane = 568.7;
double T_cr_ioctane = 543.9;
double T_b_ioctane = 543.9;
double P_cr_isooctane;

double c6_mw = 86.178;
double T_cr_3mp = 504.43;
double T_b_3mp = 336.42;
double P_cr_3mp = 31.24e+5;


double c7_mw = 100.21;
double T_cr_c7 = 540.26;
double T_b_c7 = 371.58;
double omega_c7 = 0.351;
double P_cr_c7 = 27.36e+5;

double T_cr_ndodecane = 658.0;
double T_b_ndodecane = 489.47;
double P_cr_ndodecane = 18.2e5;
double c_ndodecane_mw = 170.34;


double c16_mw = 226.45;
double T_cr_c16 = 720.60;
double T_b_c16 = 560.01;
double omega_c16 = 0.747;
double P_cr_c16 = 14.19e+5;

double h2o_mw = 18.0;
double T_cr_h2o = 647.13;
double T_b_h2o = 373.15;
double omega_h2o = 0.3449;
double P_cr_h2o = 220.55e+5;

double poly_mw = 120000.0;

double T_boil;

//air
double air_mw = 28.967;

double heat_transfer;
