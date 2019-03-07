/******************************************************************
This is code for simulation of water in n-dodecane droplet heating
to predict times to puffing

Copyright (c) 2018 Oyuna Rybdylova - All Rights Reserved
You may use, distribute and modify this code under the terms of the MIT license
Please see LICENSE file in the repository https://github.com/or38/microexplosion

If you use the software, please cite
S.S. Sazhin, O. Rybdylova, C. Crua, M. Heikal, M.A. Ismael, Z. Nissar, A. Rashid B.A. Aziz,
A simple model for puffing/micro-explosions in water-fuel emulsion droplets,
International Journal of Heat and Mass Transfer, Volume 131, 2019, Pages 815-821,
doi.org/10.1016/j.ijheatmasstransfer.2018.11.065

The details of the model implemented in the code are in the paper:
S.S. Sazhin, O. Rybdylova, C. Crua, M. Heikal, M.A. Ismael, Z. Nissar, A. Rashid B.A. Aziz,
A simple model for puffing/micro-explosions in water-fuel emulsion droplets,
International Journal of Heat and Mass Transfer, Volume 131, 2019, Pages 815-821,
doi.org/10.1016/j.ijheatmasstransfer.2018.11.065

******************************************************************/



#include "header.h"

int main(){

	int ret;
	int counter;
	

	string names;
	char resultsname[40];
	FILE * fileout;

	ret = Ini();

	//double k_l = liquid_thermcond_di(T_average);
	//double c_p_liq = liquid_Heat_Capacity(T_average);

	names = run_name;
	names += "_report.txt";
	strcpy(resultsname, names.c_str());
	fileout = fopen(resultsname, "w");
	time_step = 0.e-15;
	Fo_w = 0.0e-15;
	ret = Report_main(fileout);
	ret = Report_T(time);

	counter = 0;
	int flag = 1;
	Fo_w = Fo_ts;
//	while ((Fo_w < time_end)&&flag)
	//while ((Thickness > 0.0))
	{
		time_step = time_end;
		ret = Heat_Mass();

		//if (abs(Fo_w - 0.0001) < 1.e-6) {
			ret = Report_T(time_step);
		//}
		/*if (abs(Fo_w - 0.001) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (abs(Fo_w - 0.01) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (abs(Fo_w - 0.1) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (abs(Fo_w - 1.) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (abs(Fo_w - 10.0) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (abs(Fo_w - 100.) < 1.e-6) {
			ret = Report_T(time_step);
		}
		if (Temperature_water[N_thresh]>T_b_h2o)
			flag = 0;*/
		

//		ret = Report_T(time_step);
		ret = Report_main(fileout);
		ret = New_thickness();

	}
	//ret = Report_main(fileout);

	fclose(fileout);
	return 0;
}
int Ini(){
	string names;
	FILE * inputfile;
	char somename[50];
	int scan, i;

	cout << "Enter the name of file with input data" << endl;
	cin >> run_name;
	names = run_name;
	names += ".txt";
	strcpy(somename, names.c_str());

	if (!(inputfile = fopen(somename, "r"))) return 1;

	scan = fscanf(inputfile, "%lf", &T_gas);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &P_gas);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &Thickness_ini);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &T_ini);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &water_vol_ini);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &Fo_ts);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%lf", &time_end);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%d", &N_Lambda);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%d", &N_INT);
	fgets(somename, 50, inputfile);
	scan = fscanf(inputfile, "%d", &iprint);
	fgets(somename, 50, inputfile);
	fclose(inputfile);

	N_thresh = N_INT;// -N_INT / 20;
	//micrometres to meters
	Thickness_ini /= 1.e6;
	Thickness = Thickness_ini;

	//water radius
	water_Rad = pow(water_vol_ini, 1.0 / 3.0)*Thickness_ini;
	rad0 = water_Rad;

	Temperature_water.resize(N_INT + 1, 0.0e-15);
	Temperature_fuel.resize(N_INT + 1, 0.0e-15);
	Temperature_water1.resize(N_INT + 1, 0.0e-15);
	Temperature_fuel1.resize(N_INT + 1, 0.0e-15);
	//Initial film temperature is set to be the same as the wall temperature
	Delta_Rw = water_Rad / ((double)(N_INT));
	Delta_Rf = (Thickness - water_Rad) / ((double)(N_INT));

	for (i = 0; i < N_INT + 1; i++) Temperature_water[i] = T_ini;
//		Temperature_water[i] = T_ini + 1.0 / water_Rad * ((double) i)*Delta_Rw;
	for (i = 0; i < N_INT + 1; i++) Temperature_fuel[i] = T_ini;
//		Temperature_fuel[i] = T_ini + 1.0 - water_Rad / (Thickness - water_Rad) + 1.0 / (Thickness - water_Rad) * (water_Rad + ((double)i)*Delta_Rf);

	T_w_av = T_ini;
	T_f_av = T_ini;

	//initial mass and vol fractions
	scan = Update_vol_fracs();
	scan = Update_mass_fracs();
	
	T_eff = T_ini;

	lambda.resize(N_Lambda, 0.0e-15);
	series.resize(N_Lambda, 0.0e-15);

	Delta_Rw = water_Rad / ((double)(N_INT));
	Delta_Rf = (Thickness - water_Rad) / ((double)(N_INT));
	
	Nu = 2.0;
	T_ref = (T_gas + 2.0 * T_ini) / 3.0;

	double cw = heat_capacity_water(T_w_av);
	double rw = density_water(T_w_av);
	double kw = thermal_conductivity_water(T_w_av);
	kappa0 = kw / cw / rw;
	time_step = Fo_w*rad0*rad0 / kappa0;

	return 0;
}
int find_Lambda()
{
	//FILE * fout;
	int i, indexi, flag;
	const int Max_iterations = 200000;
	double lambda_left, lambda_right, f_left, f_right, lambda_mid, f_mid, count1, count2, aperiod, period2, saved_point, period_left, period_right;
	double kw, aw, kf, af, period;
	double conv_crit = 1.e-14;
	double step = 1.e-7;

	for (i = 0; i < N_Lambda; i++) lambda[i] = -1.0;

	kw = thermal_conductivity_water(T_w_av);
	kf = thermal_conductivity_diesel(T_f_av);
	aw = sqrt(density_water(T_w_av)*heat_capacity_water(T_w_av) / thermal_conductivity_water(T_w_av));
	af = sqrt(density_diesel(T_f_av)*heat_capacity_diesel(T_f_av) / thermal_conductivity_diesel(T_f_av));

	count1 = 0;
	count2 = 0;
	flag = 0;

	period = fmin(PI/(aw*water_Rad), PI/(af*(Thickness-water_Rad)));
	period2 = fmax(PI / (aw*water_Rad), PI / (af*(Thickness - water_Rad)));
	lambda_left = step;
	lambda_right = step;
	saved_point = 0.0 - step;
	for (i = 0; i < N_Lambda;)
	{
		indexi = 0;
		lambda_left = saved_point + 2.0*step;
		lambda_right = fmin(((double)(count1 + 1))*period - step, ((double)(count2) + 1)*period2 - step);
		saved_point = lambda_right;
		if (fabs(lambda_right - ((double)(count2)+1)*period2 + step) < 1.1*step)
			count2++;
		else count1++;

		f_left = aw*kw / tan(lambda_left*aw*water_Rad) + af*kf / tan(lambda_left*af*(Thickness - water_Rad)) - (kw - kf) / water_Rad / lambda_left;
		f_right = aw*kw / tan(lambda_right*aw*water_Rad) + af*kf / tan(lambda_right*af*(Thickness - water_Rad)) - (kw - kf) / water_Rad / lambda_right;
		if (f_left*f_right < 0.0)
		{
			while ((lambda_right - lambda_left > conv_crit) && (indexi<Max_iterations))
			{
				lambda_mid = (lambda_left + lambda_right)*0.5;
				f_mid = aw*kw / tan(lambda_mid*aw*water_Rad) + af*kf / tan(lambda_mid*af*(Thickness - water_Rad)) - (kw - kf) / water_Rad / lambda_mid;
				if (f_left*f_mid < 0.0)
				{
					lambda_right = lambda_mid;
					//f_right = aw*kw / tan(lambda_right*aw*water_Rad) + af*kf / tan(lambda_right*af*(Thickness - water_Rad)) - (kw - kf) / water_Rad / lambda_right;
					f_right = f_mid;
				}
				else
				{
					lambda_left = lambda_mid;
					f_left = f_mid;
					//f_left = aw*kw / tan(lambda_left*aw*water_Rad) + af*kf / tan(lambda_left*af*(Thickness - water_Rad)) - (kw - kf) / water_Rad / lambda_left;
				}
				indexi++;
			}
			lambda[i] = lambda_left;
		}
		i++;
	}
	//Message("\n");
	//	Message("Lambdas are calculated\n");
	//	fout = fopen("lambda.txt", "w");
	//	for (i = 0; i < N_Lambda; i++)
	//		fprintf(fout, "%d\t%20.19f\t%e\n", i + 1, lambda[i], lambda[i] * cos(lambda[i]) + h_0*sin(lambda[i]));
	//	fclose(fout);
	//	Message("Lambdas are printed in lambda.txt\n");
	return 0;
}
int Heat_Mass()
{
	// surface molar and mass fractions calculations for n-dodecane
	double xs_tot = 0.e-15;
	double P_sat = 0.0; // Saturation pressure
	double T_s = Temperature_fuel[N_INT];
	T_s = T_gas;
	int i, j;

	//Spalding head and mass trasnfer numbers
	//-------------------------------------------------------------------------
	// Calculate total evaporation rate
	T_ref = (T_gas + 2.0 * Temperature_fuel[N_INT]) / 3.0; //Sazhin, Progress in Energy and Combustion Science 32 (2006) 162–214 
	double rho_gas_s = P_gas / (R_air*T_ref); // ideal gas law
	//rho_gas_s = 1.225;
//	double c_p_die = (0.2979 + 1.4394*(T_ref / 300.0) - 0.1351*(T_ref / 300.0)*(T_ref / 300.0))*1000.0; //n-dodecane vapour Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
	double visc = Viscosity_gas(T_ref);
	double c_p_vap = vapour_heatcap(T_ref);

	double kgas = T_Conductivity_gas(T_ref);
	//heat_transfer = 2000.0;//constant
	//heat_transfer = 15.0;//constant
	
	Sh_Star = 2.0;
	Nu_star = 2.0;

	Nu = Nu_star;

	double c_p_gas = Heat_Capacity_gas(T_ref);
	double Pr = visc *c_p_gas / (kgas);  //Prandtl number

	//Temperature distribution

	delta_dot = 0.0;


	T_eff = T_gas;

	for (i = 0; i < N_Lambda; i++) { lambda[i] = -1.0; }
	int err = find_Lambda();
	//err = print_lambda();
	for (i = 0; i < N_Lambda; i++) { series[i] = 0.e-15; }

	double cw, rw, aw, cf, rf, af, kw, kf;
	cw = heat_capacity_water(T_w_av);
	rw = density_water(T_w_av);
	kw = thermal_conductivity_water(T_w_av);
	aw = sqrt(cw*rw / kw);

	cf = heat_capacity_diesel(T_f_av);
	rf = density_diesel(T_f_av);
	kf = thermal_conductivity_diesel(T_f_av);
	af = sqrt(cf*rf/kf);

	double T_s1;
	double aaaa = 1.0;
	double T0;
	check_A = 0.e-15;
	check_Tc = 0.e-15;
	check_Tr = 0.e-15;
//	double tt;
	double qw_n, qf_n, v_n2;
	T_s = aaaa*T_s + (1.0 - aaaa)*T_gas;
	for (i = 0; i < N_INT + 1; i++) Temperature_water1[i] = Temperature_water[i];
	for (i = 0; i < N_INT + 1; i++) Temperature_fuel1[i] = Temperature_fuel[i];
	do{
		T_s1 = T_s;
		for (i = 0; i < N_INT + 1; i++) Temperature_water[i] = Temperature_water1[i];
		for (i = 0; i < N_INT + 1; i++) Temperature_fuel[i] = Temperature_fuel1[i];

		for (i = 0; i < N_Lambda; i++)  {
			v_n2 = 0.5*cw*rw*water_Rad / sin(lambda[i] * aw*water_Rad) / sin(lambda[i] * aw*water_Rad) + 0.5*cf*rf*(Thickness - water_Rad) / sin(lambda[i] * af*(Thickness - water_Rad)) / sin(lambda[i] * af*(Thickness - water_Rad)) - 0.5*(kw - kf) / water_Rad / lambda[i] / lambda[i];
			
			T0 = T_gas - T_ini;
			qw_n = 0.e-15;
			qw_n = T0*cw*rw / v_n2 / lambda[i] / lambda[i] / aw / aw*(lambda[i] * aw * water_Rad / tan(lambda[i] * aw * water_Rad) - 1.0);

			qf_n = 0.e-15;
			//qf_n = T0 * cf*rf / v_n2 / lambda[i] / lambda[i] / af / af*(-2.0*lambda[i] * af*water_Rad / tan(lambda[i] * af*(water_Rad - Thickness)) - 1.0 / sin(lambda[i] * af*(Thickness - water_Rad)) - 1.0);
			qf_n = T0 * cf*rf / v_n2 / lambda[i] / lambda[i] / af / af * (0.0 - lambda[i] * af*water_Rad / tan(lambda[i] * af*(water_Rad - Thickness)) + lambda[i] * af*Thickness / sin(lambda[i] * af*(Thickness - water_Rad)) + 1.0);
			

			//T_w_av = Temperature_water[0];
			//tt = (T_s - Temperature_water[N_INT])*cw*rw / lambda[i] / aw*(water_Rad / tan(lambda[i] * aw*water_Rad) - 1.0 / lambda[i] / aw);
			series[i] = 0.0-(qw_n + qf_n)*exp(0.0 - lambda[i] * lambda[i] * time_step);// / v_n2;
			//series[i] = exp(0.0 - lambda[i] * lambda[i] * time_step) / v_n2 / lambda[i] * (T_s - T_w_av) * sqrt(kw * cw * rw) * (water_Rad / tan(lambda[i] * aw * water_Rad) - 1.0 / lambda[i] / aw);

			//check_A = (water_Rad / tan(lambda[i] * aw*water_Rad) - 1.0 / lambda[i] / aw) / v_n2 / tan(lambda[i] * af*(water_Rad - Thickness));
			//check_Tc += check_A * exp(0.0 - lambda[i] * lambda[i] * time_step);
			//check_Tr += 0.0;
		}

		//check_Tc = (T_s - T_w_av) * kw * aw * af * check_Tc;
		//check_chi = 1.0 + check_Tc / (T_s - Temperature_water[N_INT]);
		
		for (j = 0; j < N_INT + 1; j++) { Temperature_water[j] = T_s; }
		for (i = 0; i < N_Lambda; i++)  {
			Temperature_water[0] += series[i] * lambda[i] * aw / sin(lambda[i] * aw * water_Rad);
			for (j = 1; j < N_INT + 1; j++) Temperature_water[j] += series[i] * sin(lambda[i] * aw*((double)j)*Delta_Rw) / sin(lambda[i] * aw*water_Rad) / (((double)j)*Delta_Rw);
		}

		for (j = 0; j < N_INT + 1; j++) { Temperature_fuel[j] = T_s; }
		for (i = 0; i < N_Lambda; i++)  {
			for (j = 0; j < N_INT + 1; j++) Temperature_fuel[j] += series[i] * sin(lambda[i] * af*(Thickness - ((double)j)*Delta_Rf - water_Rad)) / sin(lambda[i] * af*(Thickness - water_Rad)) / (water_Rad + ((double)j)*Delta_Rf);
		}
		// Now we know temperature at each layer
		T_s = Update_Ts();
	} while (fabs(T_s - T_s1) > 1.e-8);

	//T_s = Temperature_fuel[N_INT];

	// Re-calculate droplet avarage temperature T_av
	T_w_av = Temperature_water[N_INT] * water_Rad * water_Rad;
	for (int j = 1; j < N_INT; j += 2) {
		T_w_av += 4.0 * Temperature_water[j] * (((double)j)*Delta_Rw)*(((double)j)*Delta_Rw);
	}
	for (int j = 2; j < N_INT; j += 2) {
		T_w_av += 2.0 * Temperature_water[j] * (((double)j)*Delta_Rw)*(((double)j)*Delta_Rw);
	}
	T_w_av = T_w_av*Delta_Rw / water_Rad / water_Rad / water_Rad;

	T_f_av = Temperature_fuel[0] * water_Rad * water_Rad + Temperature_fuel[N_INT] * Thickness * Thickness;
	for (int j = 1; j < N_INT; j += 2) {
		T_f_av += 4.0 * Temperature_fuel[j] * (water_Rad + ((double)j)*Delta_Rf)*(water_Rad + ((double)j)*Delta_Rf);
	}
	for (int j = 2; j < N_INT; j += 2) {
		T_f_av += 2.0 * Temperature_fuel[j] * (water_Rad + ((double)j)*Delta_Rf)*(water_Rad + ((double)j)*Delta_Rf);
	}
	
	T_f_av = T_f_av*Delta_Rf / (Thickness * Thickness * Thickness - water_Rad * water_Rad * water_Rad);

	boundary1 = kw*(Temperature_water[N_INT] - Temperature_water[N_INT - 1]) / Delta_Rw;
	boundary2 = kf*(Temperature_fuel[1] - Temperature_fuel[0]) / Delta_Rf;
	//boundary_cond1 = heat_transfer*(T_eff - Temperature[N_INT]);
	//boundary_cond2 = k_l * (Temperature[N_INT] - Temperature[N_INT - 1]) / Delta_R/Thickness;


/*	string names;
B	FILE * fileout;
	stringstream ss;

	names = run_name;
	names += "_T_";
	ss << time;
	names += ss.str();
	names += ".txt";

	strcpy(resultsname, names.c_str());
	fileout = fopen(resultsname, "w");

	for (int j = 0; j < N_INT + 1; j++) fprintf(fileout, "%e\t", Temperature[j]);
	fprintf(fileout, "\n");

	fclose(fileout);*/

	return 0;
}
int New_thickness()
{
	int i;
	Fo_w += Fo_ts;

	for (i = 0; i < N_INT + 1; i++) Temperature_water[i] = T_ini;
	for (i = 0; i < N_INT + 1; i++) Temperature_fuel[i] = T_ini;

	T_w_av = T_ini;
	T_f_av = T_ini;

	T_eff = T_ini;

	Nu = 2.0;
	T_ref = (T_gas + 2.0 * T_ini) / 3.0;

	return 0;
}
double Update_Ts()
{
	double ts = T_gas;//Temperature_fuel[N_INT];
	return ts;
}

//air 
//from http://www.thermopedia.com/content/553/
double Viscosity_gas(double T)
{
	//return 1.85e-5; //at 300 K
	return (4.5608 + 0.70077*T - 0.00037287*T*T + 0.000000090437*T*T*T)*1.e-6;
	//return 4.0e-5;
}
double T_Conductivity_gas(double T)
{
	//return 0.0263; // at 300 K
	return -0.00038603 + 0.00010311*T - 0.000000054199*T*T + 0.000000000017429*T*T*T;
	//return 0.07;
}
double Heat_Capacity_gas(double T)
{
	//return 1005.0;//at 300 K
	return (29.643 - 0.0051373*T + 0.000013106*T*T - 0.0000000048325*T*T*T) / air_mw*1.e3;
	//return 1127.0;
}
//vapour
//liquid
double liquid_density(double T)
{
	//average density
	double den1, den2, den;
	den1 = density_water(T);
	den2 = density_diesel(T);
	den = vol_w*density_water(T_w_av) + vol_f*density_diesel(T_f_av);
	return den;
}
double Saturation_Pressure_3mp(double T)
{
	/*double A = 50.3422;
	double B = -3.2789e3;
	double C = -1.6111e1;
	double D = 7.4260e-3;
	double E = -9.1804e-14;
	A = 35.2848;
	B = -2.6773e3;
	C = -9.8546;
	D = 2.2352e-11;
	E = 4.0277e-6;
	double p = pow(10.0, A + B / T + C*log10(T) + D*T + E*T*T);
	return pow(10.0, A + B / T + C*log10(T) + D*T + E*T*T)*101325.0 / 760.0; //Yaws 1997, as in Zhang, 2017*/
	double Tr = T / T_cr_3mp;
	double tau = 1 - Tr;
	double f0, f1, f2;
	f0 = (-5.97616*tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5.0)) / Tr;
	f1 = (-5.03365*tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5.0)) / Tr;
	f2 = (-0.64771*tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5.0)) / Tr;
	double omega1 = 0.274;
	return exp(f0 + f1 * omega1 + f2 * omega1* omega1)*1.e+5;
}
double Saturation_Pressure_ic8(double T)
{
/*	double A = 50.3422;
	double B = -3.2789e+3;
	double C = -1.6111e+1;
	double D = 7.4260e-3;
	double E = -9.1804e-14;
	double p = pow(10.0, A + B / T + C*log10(T) + D*T + E*T*T);
	return pow(10.0, A + B / T + C*log10(T) + D*T + E*T*T)*101325.0 / 760.0; //Yaws 1997, as in Zhang, 2017*/
	double Tr = T / T_cr_ioctane;
	double tau = 1 - Tr;
	double f0, f1, f2;
	f0 = (-5.97616*tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5.0)) / Tr;
	f1 = (-5.03365*tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5.0)) / Tr;
	f2 = (-0.64771*tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5.0)) / Tr;
	double omega1 = 0.303;
	return exp(f0 + f1 * omega1 + f2 * omega1* omega1)*1.e+5;

}
double Saturation_Pressure_c7(double T)
{
	double Tr = T / T_cr_c7;
	double tau = 1 - Tr;
	double f0, f1, f2;
	if (T > 0.99*T_cr_c7) { Tr = 0.99; tau = 0.01; }
	f0 = (-5.97616*tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5.0)) / Tr;
	f1 = (-5.03365*tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5.0)) / Tr;
	f2 = (-0.64771*tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5.0)) / Tr;
	return exp(f0 + f1 * omega_c7 + f2 * omega_c7* omega_c7)*1.e+5;
}
double Saturation_Pressure_c16(double T)
{
	double Tr = T / T_cr_c16;
	double tau = 1 - Tr;
	double f0, f1, f2;
	if (T > 0.99*T_cr_c16) { Tr = 0.99; tau = 0.01; }
	f0 = (-5.97616*tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5.0)) / Tr;
	f1 = (-5.03365*tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5.0)) / Tr;
	f2 = (-0.64771*tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5.0)) / Tr;
	return exp(f0 + f1 * omega_c16 + f2 * omega_c16* omega_c16)*1.e+5;
}
double Saturation_Pressure_h2o(double T)
{
	double Tr = T / T_cr_h2o;
	double tau = 1 - Tr;
	double f0, f1, f2;
	if (T > 0.99*T_cr_h2o) { Tr = 0.99; tau = 0.01; }
	f0 = (-5.97616*tau + 1.29874 * pow(tau, 1.5) - 0.60394 * pow(tau, 2.5) - 1.06841 * pow(tau, 5.0)) / Tr;
	f1 = (-5.03365*tau + 1.11505 * pow(tau, 1.5) - 5.41217 * pow(tau, 2.5) - 7.46628 * pow(tau, 5.0)) / Tr;
	f2 = (-0.64771*tau + 2.41539 * pow(tau, 1.5) - 4.26979 * pow(tau, 2.5) + 3.25259 * pow(tau, 5.0)) / Tr;
	return exp(f0 + f1 * omega_h2o + f2 * omega_h2o* omega_h2o)*P_cr_h2o;
}
double Saturation_Pressure_diesel(double T)
{
	double p;
	p = exp(8.1948 - 7.8099*(300.0 / T) - 9.0098*(300.0 / T)*(300.0 / T))*1.e5;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
	if (T > 0.99*T_cr_ndodecane) {
		p = p*exp(15.0*(T / 0.99 / T_cr_ndodecane - 1.0));
	}
	//return exp(8.1948 - 7.8099*(300.0 / T) - 9.0098*(300.0 / T)*(300.0 / T))*1.e5;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
	return p;
}

double Evaporation_Heat_ic8(double T)
{
	return 49.32456*pow(1 - T / T_cr_ioctane, 0.382229) / c8_mw * 1.e6;
}
double Evaporation_Heat_h2o(double T)
{
	return 54.0*pow(1 - T / T_cr_h2o, 0.34) / h2o_mw * 1.e6;
}
double Evaporation_Heat_3mp(double T)
{
	return 1.093*R_UNIVERSAL*T_b_3mp*(log(P_cr_3mp*1.e-5) - 1.013) / (0.93 - T_b_3mp / T_cr_3mp);
}
double Evaporation_Heat_c7(double T)
{
	if (T>0.99 * T_cr_c7) return 49.73*pow(1 - 0.99, 0.386) / c7_mw * 1.e6;
	return 49.73*pow(1 - T / T_cr_c7, 0.386) / c7_mw * 1.e6;
}
double Evaporation_Heat_c16(double T)
{
	if (T>0.99 * T_cr_c16) return 96.68*pow(1 - 0.99, 0.4220) / c16_mw * 1.e6;
	return 96.68*pow(1 - T / T_cr_c16, 0.4220) / c16_mw * 1.e6;
}
double Evaporation_Heat_ndodecane(double T)
{
	double l;
	l = 37.44*pow(T_cr_ndodecane - T, 0.38)*1000.0;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
	if (T>0.99*T_cr_ndodecane) {
		l = 37.44*pow(5.0, 0.38)*1000.0;
	}
	return l;

}

double Diffusion_coef_Vap(double T, double p)
{
	return 0.527*pow(T / 300.0, 1.583) / p; //n-dodecane
}
double vapour_heatcap(double T)
{
	return (0.2979 + 1.4394*(T / 300.0) - 0.1351*(T / 300.0)*(T / 300.0))*1000.0;
}

double density_water(double T)
{
//	return 600;
	double den;
	den = 0.325 * pow(0.27, -pow(1.0 - (T / T_cr_h2o), 0.23))*1000.0;
	return den;
//Yaws
}
double heat_capacity_water(double T)
{
//	return 2830;
	double TC1;
	TC1 = -22.417 + 8.7697e-1 * T - 2.5704e-3 * T * T + 2.4838e-6 * T * T * T;
	TC1 = TC1 / h2o_mw * 1000.0;
	return TC1;
	//Yaws
}
double thermal_conductivity_water(double T)
{
//	return 0.145;
	double TC1;
	TC1 = -3.5667e-1 + 5.0570e-03*T - 6.1071e-6*T*T;
	return TC1;//Yaws
}

double density_diesel(double T)
{
//	return 23.8;
	return 744.11 - 0.771*(T - 300.0);
}
double heat_capacity_diesel(double T)
{
//	return 1120.0;
	return (2.18 + 0.0041*(T - 300.0))*1000.0;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
}
double thermal_conductivity_diesel(double T)
{
	//return 6.1;
	return 0.1405 - 0.00022*(T - 300.0); //Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
}

int Update_vol_fracs()
{
	vol_w = pow(water_Rad / Thickness, 3.0);
	vol_f = 1.0 - vol_w;
	return 0;
}
int Update_mass_fracs()
{
	Y_w = vol_w*density_water(T_w_av) / (vol_w*density_water(T_w_av) + vol_f*density_diesel(T_f_av));
	Y_f = 1.0 - Y_w;
	return 0;
}

//Technical
int Report_T(double ttime)
{
	string names;
	char resultsname[40];
	FILE * fileout;
	stringstream ss;

	names = run_name;
	names += "_T_";
	ss << ttime;
	names += ss.str();
	names += ".txt";

	strcpy(resultsname, names.c_str());
	fileout = fopen(resultsname, "w");

	for (int i = 0; i < N_INT + 1; i++) fprintf(fileout, "%e\t%e\n", (((double)i)*Delta_Rw) / water_Rad, Temperature_water[i]);
	for (int i = 1; i < N_INT + 1; i++) fprintf(fileout, "%e\t%e\n", (water_Rad + ((double)i)*Delta_Rf) / water_Rad, Temperature_fuel[i]);


	fclose(fileout);
	return 0;
}
int Report_main(FILE * ff)
{
	double kg = 0.061 / 1120.0 / 23.8;
	fprintf(ff, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", time_step, Fo_w, Thickness, T_w_av, Temperature_water[N_thresh], (boundary1 - boundary2)*2.0 / (boundary1 + boundary2), boundary1, boundary2, check_chi);
	return 0;
}
int print_lambda()
{
	double cw, rw, aw, kw;
	cw = heat_capacity_water(T_w_av);
	rw = density_water(T_w_av);
	kw = thermal_conductivity_water(T_w_av);
	aw = sqrt(cw*rw / kw);
	FILE *file;
	file = fopen("lambda.txt", "w");
	for (int i = 0; i < N_Lambda; i++) fprintf(file, "%f\n", lambda[i] * aw * water_Rad/PI);
	fclose(file);
	return 0;
}
