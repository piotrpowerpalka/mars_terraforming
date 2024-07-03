/**
* Name: MarsEnergyModel
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
* 

*/


model MarsEnergySoilBucket


/* Insert your model definition here */

global {
	float f <- 10e-4; 			// [sec-1] - Coriolis parameter
	float R <- 0.19e-3; 		// [J kg-1 K-1] gas constant = 0.19×107 ergs gm-1 K-1
	float S <- 2e-7;			// [K *m-1] static stability
	float L <- 5.3e6;			// [m] equator-pole distance
	float ro <- 2.4733;		// [kg m-3] CO2 density in 220
	float south_dsg <- 0.0; 
	float north_dsg <- 0.0; 
	float total_sg  <- 0.0;  // total sun glasses
	float solar_const_modifier <- 1.0;
	float k_reg <- 1.0; // 0.0837;  // thermal conductivity of soil [W m-1 K-1]
	float deepSoil <- 3.47; // 60.0; // [m] - grubość gruntu głębokiego
	float initDeepGrdTemp <- 212.0; // [K] - initial deep ground temp.
	float kreg_mod <- 1.0;		// k regolith modifier for areas near poles
	float kreg_lat <- 90.0;		// from what lattitude, the kreg_mod will be modified
	float ndsg_lat <- 90.0;
	float sdsg_lat <- -90.0;
	float spCO2defreze <- 12.5e-6;	// south pole CO2 deffreze param
	
	float GHG_A <- 0.004;
	float GHG_B <- 0.4551;
	
	float totalCO2 <- 0.007; // [bar] total amount of CO2 on Mars
	float sumFrozen <- 0.0;
	float atmCO2 <- totalCO2 - sumFrozen;
	//float sumFrozen update: cell sum_of(each.frozenCO2); 
	//float atmCO2 update: totalCO2 - sumFrozen;  // [bar] mean CO2 pressure
	
	int aspect_mode <- 1;				// tryb wyświetlania
	int sollon_number <- 0;
    int model_number <- 3;
    int selected_cell <- 1788;
    list<int> selected_cells <- [1788, 2859];
    int hadleyModel <- 0; // 0: nasz model - prosty, 1: model Alison
    int deepGroundModel <- 0; // 0: brak, 1: jest
    
    int numOfHexes <- 4002;
    float Tincrease <- 0.0;
    
    int GHGregress <- 0;
    //float Gamma_HT <- 10/1000; // Kelvin / m - potential temperature parameter
    float Gamma_HT <- ga/cCO2; // 0.004751 wzór wg. Handbook on Atmospheric Diffusion
     
	float ecc <- 0.09341233; // eccentricity 
	
    float tempFrozen <- 0.0;
    float tempTemp <- 0.0;
    float tempCO2 <- 0.0;
    float tempEnergy <- 0.0;
    float eddyTemp <- 0.0;
    float hadTemp <- 0.0;
    float latentEn <- 0.0;
	float greenEn <- 0.0;
	float radEn <- 0.0;
	float insolEn <- 0.0;
	float prevEn <- 0.0;
	float eddyEn <- 0.0;
	float hadEn <- 0.0;
	float tempSat <- 0.0;
	float tempTsat <- 0.0;
	float grEn <- 0.0;
	float sgrTs <- 0.0;
	
	float fitness <- 1e12;
	float yty_changes_RMSE <- 1e12;
	float yty_diff <- 1e12;
	
	float mape_tv1 <- 0.0;
	float mape_pv1 <- 0.0;
	float mape_tv2 <- 0.0;
	float mape_pv2 <- 0.0;
	float rmse_tv1 <- 0.0;
	float rmse_pv1 <- 0.0;
	float rmse_tv2 <- 0.0;
	float rmse_pv2 <- 0.0;
	
	float rmse_avgt <- 0.0;
	float mape_avgt <- 0.0;
	
	float rmse_avgp <- 0.0;
	float mape_avgp <- 0.0;
	
	
	float tauDust <- 0.0;
	float albedoGlobal <- 0.2;	
	float albedoIce <- 0.5;
		
	string outdir <- "../results/";
	    
	int log_every_sol <- 1;
	bool log_output <- false;
	bool log_csv <- false;
	bool log_vik <- false;
	bool log_calib <- false;
	int log_first_year <- 0;   // first year that is logged to output
	bool sunGlass_Ice <- false;		// on constant set of hexes: 0, put sunglasses when ice only: 1
	int sun_glass_param <- 0; // 0: czy "okulary" sa ustawione: 0: na stałych hexach (numery), 1: zależnie od kreg_lat, 2: czy zależą od pokrywy lodowej (false)
	int sun_glass_mode <- 0;  // 0: okulary takie same w zadanym rejonie, 1: sinus od 0 (kreg_lat) do 1 (biegun)
	bool variableTauDust <- false;
	float scaleHeight <- 11600; //[m] 
	float constrTemp <- 50.0;		// [K] temperature constr.
	bool albedoFile <- false;
	
	float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
    float ga <- 3.72076; 			// Mars gravitional acceleration [ms-2] [Nkg-1]
	float AU <- 149597870700.0;		// astronomical unit [m]
	float sunTemperature <- 5780.0; 	// sun temperature [K]
	float sunRadius <- 695700000.0;	// sun radius [m]
	float marsRadius <- 3396200.0;	// mars radius [m]
	float marsArea <- 1.448e14;		// mars area [m2]\
	float S0 <- 589;		// solar constant 
	float e <- 0.09341233;		// mimośród (spłaszczenie orbity)
	float nachylenieOsi <- 24.936;

	float cCO2 <- 783.0; // specific heat cappacity for CO2 [J kg-1 K-1] - in 220 [K] https://en.wikipedia.org/wiki/Carbon_dioxide
	float cSoil <- 980.0;  // specific heat cappacity for soil (desert sand) [J kg-1 K-1] https://aip.scitation.org/doi/pdf/10.1063/1.4949109
	float soilDepth <- 2.0; // [m]
	float soilDensity <- 1600; // [kg m-3]
	float emissivity <- 0.9; // [W m-2]
	
	float delta_t <- 88775.0;		// 1 martian sol [s]
	float delta_h2 <-  marsArea / numOfHexes;	// area of single hex
	float delta_h <- sqrt(delta_h2 / #pi);			// radius of single hex (approx.)
	float Kcoeff <- 10000.0;
	float K_ <- Kcoeff * delta_t / delta_h2 / step_sol; // K * delta_t / delta_h2
	float Wind <- 0.0;
	float initFrozen <- 0.0;		// część początkowa zamrożonego CO2
		
   	float hadleyN <- 0.0;
   	float hadleyS <- 0.0;
   	//float hadleyParam <- 5e-9;
   	float hadleyParam <- 0.000000005;
   		
	int martianYear <- 668; // martian year [sol]
	float Pa2bar <- 1e5;			// pascal to bar
	
	// albedo
	file csv_albedo <- csv_file("../includes/heksy_albedo.csv",";",float,true);
	matrix<float> mat_albedo <- matrix<float>(csv_albedo.contents);
	
	file csv_hadley <- csv_file("../includes/hadley_heat_transport.csv",",",float,true);
	matrix<float> mat_hadley <- matrix<float>(csv_hadley.contents);
	
	file csv_v1 <- csv_file("../includes/v1.csv",";",float,true);
	matrix<float> vik1 <- matrix<float>(csv_v1.contents);
	
	file csv_v2 <- csv_file("../includes/v2.csv",";",float,true);
	matrix<float> vik2 <- matrix<float>(csv_v2.contents);
	
	file csv_calib <- csv_file("../includes/calib_data.csv",";",float,true);
	matrix<float> calib <- matrix<float>(csv_calib.contents);
	
		
	
	int step_sol <- 1;	// number of iterations per sol
	int sol_year <- 0; // update: (cycle/(2 * step_sol)) mod martianYear; // number of sol in martian year
    int sol_lon  <- 0; // update: int(sol_year / 668.0 * 360.0); // recalculation sol no -> solar logitude
    int year <- 0; // update: int( cycle/(2 * step_sol)) - sol_year)/martianYear;   
    
	float mean_greenhouse_temp 	<- 0.0; // update: cell mean_of (each.Ts);  // mean greenhouse temperature
	float var_greenhouse_temp 	<- 0.0;
	float mean_equatorial_temp 	<- 0.0; //update: equ_list mean_of (each.Ts);
	float mean_Npole_temp 		<- 0.0;// update: N_poles_list mean_of(each.Ts); 
	float mean_Spole_temp 		<- 0.0; //update: S_poles_list mean_of(each.Ts); 
	
	int biol_habitable_hex <- 0;		// number of hexes biologically habitable (T > -25 C)
	int biosphere_hex <- 0; 			// number of hexes unfreezed water (T > 0 C)
	
    
	float CO2freezingCoeff <- 0.00001 / step_sol; // percentage part of CO2 being sublimated / resublimated per sol
    
    list<cell> N_poles_list <- nil;
    list<cell> S_poles_list <- nil;
    list<cell> _45p_poles_list <- nil;
    list<cell> _45m_poles_list <- nil;
    list<cell> equ_list <- nil;
 	
 	matrix had_al <- 0.0 as_matrix({4002,668});
	
    init {
    	
		if (hadleyModel = 0) {
			csv_hadley <- csv_file("../includes/hadley_heat_transport.csv",",",float,true);
			mat_hadley <- matrix<float>(csv_hadley.contents);			
		}
		if (hadleyModel = 2) {
			csv_hadley <- csv_file("../includes/hudley_new.csv",",",float,true);
			mat_hadley <- matrix<float>(csv_hadley.contents);
		}
		if (hadleyModel = 1) {
			file hadley_al_csv <- csv_file("../includes/hudley_al.csv",";",float,true);
			int cnt <- 0;
			int id_hex <- -1;
			int id_sol <- -1;
			
			loop el over: hadley_al_csv {
				if (cnt = 0){
					id_hex <- int(el) - 1;		
				} if (cnt = 1){
					id_sol <- int(el);
				} if (cnt = 2){
					had_al[id_hex, id_sol] <- float(el);
				}
				cnt <- mod(cnt + 1, 3);
			}
		}


		create cell from: csv_file("../includes/fixed/out_grid_4002_0h_" + sollon_number + "_sollon.csv", ";", true) with:
    	[   id_cell::int(read("id")),
    		x::float(read("x")),
    		y::float(read("y")),
    		z::float(read("z")),     	
    		longitude::float(read("longitude")),
    		latitude::float(read("lattitude")),
    		spec::int(read("spec")),
		    n1::int(read("n1")),
    		n2::int(read("n2")),
    		n3::int(read("n3")),
    		n4::int(read("n4")),
	   		n5::int(read("n5")),
    		n6::int(read("n6")),
    		a1::int(read("a1")),
    		a2::int(read("a2")),
    		a3::int(read("a3")),
    		a4::int(read("a4")),
    		a5::int(read("a5")),
    		a6::int(read("a6")),
    		height::float(read("ex2"))
    		] {
    			next_step_co2_column <- co2_column;
    			int ntmp;
				int atmp; 
				albedo <- mat_albedo[id_cell-1];
				if (hadleyModel = 0 or hadleyModel = 2){
					hadley_ht <- mat_hadley[id_cell-1];	
					write mat_hadley[id_cell - 1];
				}
				frozenCO2 <- totalCO2 * initFrozen / numOfHexes;
								
  				location <- {longitude, latitude};
    			neigh <- [n1, n2, n3, n4, n5, n6];
				azim  <- [a1, a2, a3, a4, a5, a6];
				Ts <- temp;
				ground_delta_En <- 0.0;
				subground_Ts <- initDeepGrdTemp;
				
				if (spec = 1) {
					remove index: 5 from: neigh;
					remove index: 5 from: azim;
				}
				
				if (sun_glass_param = 0){
					if (id_cell in [1751, 1752, 1753, 1769, 1770, 1771, 1772, 1787, 1788, 1789, 1790, 1791, 1835, 1844, 1845, 1854, 1855, 1865, 1866]){
						dustSunGlasses <- south_dsg;
					} 
					if (id_cell in [2841, 2849, 49, 2858, 59, 2850, 2859, 58, 2869, 69, 2860, 2870, 68, 2871, 80, 2881, 79, 2882, 2883]) {
						dustSunGlasses <- north_dsg;
					}
				} else if (sun_glass_param = 1) {
					if (latitude < sdsg_lat) { // sdsg_lat jest ujemne
						if (sun_glass_mode = 0){
							dustSunGlasses <- south_dsg;	
						}
						if (sun_glass_mode = 1) {
							dustSunGlasses <- sin(90*(abs(latitude) + sdsg_lat)/(90.0 + sdsg_lat)) * south_dsg;
						}
						
					} 
					if (latitude > ndsg_lat){
						if (sun_glass_mode = 0){
							dustSunGlasses <- north_dsg;	
						}
						if (sun_glass_mode = 1) {
							dustSunGlasses <- sin(90*(abs(latitude) - ndsg_lat)/(90.0 - ndsg_lat)) * north_dsg;
						}
					}
				} else if (sun_glass_param = 2) {
					
					if (latitude > 0) {
						dustSunGlasses <- north_dsg;
					} else if (latitude < 0) {
						dustSunGlasses <- south_dsg;
					}
				}
				
				
				loop idx from: 0 to: length(neigh)-1 { 
					loop jdx from: 0 to: length(neigh)-2 {
						if (azim[idx] < azim[jdx]) {
							atmp  <- azim[idx];
							azim[idx] <- azim[jdx];
							azim[jdx] <- atmp;
							
							ntmp  <- neigh[idx];
							neigh[idx] <- neigh[jdx];
							neigh[jdx] <- ntmp;
						}
					}
				}
    		}
    		
    	N_poles_list <- cell where (each.latitude >= 90);
		S_poles_list <- cell where (each.latitude <= -90);
		_45p_poles_list <- cell where (each.latitude >= 44 and each.latitude <= 46);
		_45m_poles_list <- cell where (each.latitude >= -46 and each.latitude <= -44);
		
		
		equ_list   <- cell where (each.latitude >= -2 and each.latitude <= 2);

    }
    
    reflex update {
    	sumFrozen <- cell sum_of(each.frozenCO2); 
		atmCO2 <- totalCO2 - sumFrozen;  // [bar] mean CO2 pressure
		sol_year <- (cycle/(2 * step_sol)) mod martianYear; // number of sol in martian year
   		sol_lon  <- int(sol_year / 668.0 * 360.0); // recalculation sol no -> solar logitude
   		year <-  int( (cycle/(2 * step_sol) - sol_year)/martianYear );
    
    	mean_greenhouse_temp 	<- cell mean_of (each.Ts);  // mean greenhouse temperature
		mean_equatorial_temp 	<- equ_list mean_of (each.Ts);
		mean_Npole_temp 		<- N_poles_list mean_of(each.Ts); 
		mean_Spole_temp 		<- S_poles_list mean_of(each.Ts); 
    }
    
    /**
     * save data every year
     */
    reflex saving_output when: log_output and year >= log_first_year and (
    	cycle mod (2 * step_sol * log_every_sol) = 1 // or cycle mod 24 = 1
    	//cycle mod 56 = 1
    	) {    	
    	loop n over: cell   {
    		int sol <-  (cycle -1) / (2 * step_sol);
 	 			write "" + ( (cycle -1) / (2 * step_sol) ) + ";" + n.id_cell+ ";" + 
	 						n.Ts + ";" + n.pCO2 + ";"  +
    	          	 		n.frozenCO2  + ";" + n.heat_flux + ";" + 
    	          	 		(sol mod 668) + ";" + (sol - (sol mod 668))/668 + ";" +
    	          	 		n.insol + ";" + n.Tps + ";"  +
    	          	 		n.kat + ";" + n.latitude +  ";" + n.rad_en + ";" + n.green_en + ";" + n.eddy_Ts  
    	          	 		+ ";" + n.hadley_Ts + ";" + n.latentCO2 + ";" + n.ground_delta_En + ";" + n.subground_Ts
    	          	 		;  
    	    }
    	
    }
    reflex saving_vik when: log_vik and (
    	cycle mod (2 * step_sol * log_every_sol) = 1 // or cycle mod 24 = 1
    	//cycle mod 56 = 1
    	) {    	
    	loop n over: cell   {
    		int sol <-  (cycle -1) / (2 * step_sol);
    		
    		if (n.id_cell = 333 or n.id_cell = 3090){
 	 			write "" + sol + ";"  + n.id_cell+ ";" + 
	 						n.Ts + ";" + n.pCO2 + ";"  +
    	          	 		n.frozenCO2 + ";" + n.heat_flux  + ";" + 
    	          	 		(sol mod 668) + ";" + (sol - (sol mod 668))/668 + ";" +
    	          	 		n.insol+ ";" + n.Tps
    	          	 		;     
    	    }
    	}
    }
			 
			 

    reflex saving_csv when: log_csv and (
    	cycle mod (2 * step_sol * log_every_sol) = 1 // or cycle mod 24 = 1
    	) {    	
    	loop n over: cell   {
    		if (length(selected_cells) = 0 or n.id_cell in selected_cells) {
	    	    save [ 
	    	    	n.id_cell, n.co2_column, 
	    	    	n.Ts, n.insol, n.energy, n.pCO2, n.frozenCO2, n.heat_flux, n.Tps
	    	    ] to: outdir + "/year_" + int( cycle / 1336 )  + "_sol_" + (cycle / 2) mod 668 + "_cell_" + n.id_cell + ".csv" rewrite: true type: "csv";
    	    }    	   
     	}
    }
    
    
    reflex hadley_heat_transport when: cycle > 0 and cycle mod 2 = 1 and (hadleyModel = 0 or hadleyModel = 2){
    	hadleyN <- 0.0;
    	hadleyS <- 0.0;
    	loop n over: cell {
    		ask n {n.hadley_Ts <- 0.0; }
    	}
    	
    	
    	// first - get the heat from -30 -- + 30cells
    	loop n over: cell {
    		if (n.hadley_ht > 0.0 and n.latitude > 0) {
    			ask n {
    				hadley_Ts <- hadley_ht * Ts * hadleyParam;
    			}
    			hadleyN <- hadleyN + n.hadley_Ts;
    		}
    		if (n.hadley_ht > 0.0 and n.latitude < 0) {
    			ask n {
    				hadley_Ts <- hadley_ht * Ts * hadleyParam;
    			}
    			hadleyS <- hadleyS + n.hadley_Ts;
    		}
    		if (n.hadley_ht > 0.0 and n.latitude = 0) {
    			ask n {
    				hadley_Ts <- hadley_ht * Ts * hadleyParam;
    			}
    			hadleyN <- hadleyN + 0.5 * n.hadley_Ts;
    			hadleyS <- hadleyS + 0.5 * n.hadley_Ts;
    		}
    	}
    	// after - distribute heat to -90 -- -30 and +30 -- +90
    	loop n over: cell {
    		if (n.hadley_ht < 0.0 and n.latitude > 0) {
    			ask n {
    				n.hadley_Ts <- - hadley_ht * hadleyN;
    			}
    		}
			if (n.hadley_ht < 0.0 and n.latitude < 0) {
    			ask n {
    				n.hadley_Ts <- - hadley_ht * hadleyS;
    			}
    		}   		
    	}
		loop n over: cell {
    		ask n {n.hadley_Ts <- n.hadley_Ts * mean_greenhouse_temp^0.5 * (mean_equatorial_temp - (mean_Npole_temp+mean_Spole_temp)/2)^2; }
    	}
	}
	
	reflex rmse when: log_calib and cycle mod (2 * step_sol * log_every_sol) = 1 {
		float avgp <- cell mean_of(each.pCO2);
		
		rmse_avgt <- rmse_avgt + (calib[1, sol_year] - mean_greenhouse_temp)^2;
		mape_avgt <- mape_avgt + abs(calib[1, sol_year] - mean_greenhouse_temp)/calib[1, sol_year];
		
		rmse_avgp <- rmse_avgp + (calib[4, sol_year]/1e5 - avgp)^2;
		mape_avgp <- mape_avgp + abs(calib[4, sol_year]/1e5 - avgp)/(calib[4, sol_year]/1e5);
		
	}
    
    reflex reset_diff when: log_calib and sol_year = 667 and cycle mod (2 * step_sol * log_every_sol) = 1{
		mape_tv1 <- mape_tv1 / 358;
		mape_pv1 <- mape_pv1 / 358;
		
		mape_tv2 <- mape_tv2 / 367;
		mape_pv2 <- mape_pv2 / 367;
		
		rmse_tv1 <- rmse_tv1 / 358;
		rmse_pv1 <- rmse_pv1 / 358;
		
		rmse_tv2 <- rmse_tv2 / 367;
		rmse_pv2 <- rmse_pv2 / 367;
		
		rmse_avgt <- rmse_avgt / 668;
		mape_avgt <- mape_avgt / 668;
		rmse_avgp <- rmse_avgp / 668;
		mape_avgp <- mape_avgp / 668;
		
		
		
		//fitness <- mape_tv1 + mape_tv2 + mape_pv1 + mape_pv2;
		//fitness <- rmse_tv1 + rmse_tv2 + rmse_pv1 * 30000 + rmse_pv2* 30000 + rmse_avgt + rmse_avgp* 30000;
		fitness <- rmse_avgt/200 + rmse_avgp/0.0065;
		
		
		//fitness <- max(mape_tv1, mape_pv1);  
		
		if (log_calib){
			write "" + fitness + ";" + rmse_tv1 + ";" + rmse_tv2 + ";" + rmse_pv1 + ";" + rmse_pv2 + ";" + rmse_avgt + ";" + rmse_avgp;
			write "" + hadleyParam + ";" + soilDepth + ";" + totalCO2 + ";" + CO2freezingCoeff ; 
		}
	}
    
    // check whether the difference in results in small enough
   	reflex check_mtx when: log_calib and year > 1 and sol_year = 0 {
   		yty_diff <- 0.0;
		loop n over: cell {
			yty_diff <- yty_diff + n.mtx_diff/4002;
		}
	}
}



species cell parallel: true {
	list<cell> neighbours;
	
	int id_cell;
	float x;
	float y;
	float z;
	int n1;
	int n2;
	int n3;
	int n4;
	int n5;
	int n6;
	int a1;
	int a2;
	int a3;
	int a4;
	int a5;
	int a6;
	float temp <- 140.0;
	float energy <- 0.0; // energia w równaniu energy balance
	float prev_energy <- 0.0;
	
	matrix temp_mtx <- 0.0 as_matrix({2,668});  // to compare the results with the previous year
	matrix pres_mtx <- 0.0 as_matrix({2,668});  // to compare the results with the previous year
	float mtx_diff <- 1e12;
	
	float longitude;
	float latitude;
	float height;
	float zonal_wind;
	float merid_wind;
	float dustSunGlasses <- 0.0;

	float atm_press;
	int spec;
	
	float co2_column;
	float next_step_co2_column;
 	float next_step_Ts;
 	float neigh_Ts_sum;
 	
 	float div_Ts <- 0.0;
 	float heat_flux <- 0.0;
 	
	float co2_atmpres_share;
	float sunMarsDist <- 1.555252;
	float distanceFromPlanetCenter;
	float frozenCO2;
    
	float S <- 1.0; //Insolation factor 
	float albedo; // albedo
	float solarZenithAngle; // solar zenith angle
	
	
	// meanCO2 - frozenTotal 
	// pTotal = pA + froz + poles
	float pCO2 <- 0.0; // update: atmCO2 * exp(- height / scaleHeight); // 0.006 * exp(- height / scaleHeight); // Pressure in [bar]
	
	float PsatCO2 <- 0.0; // update: 1.2264e7 * exp(-3167.8 / Ts); //[bar] za Reference this from Eq 19 in Fanale et al. (1982) Fanale, F.P., Salvail, J.R., Banerdt, W.B. and 
													// Saunders, R.S., 1982. Mars: The regolith-atmosphere-cap system and climate 
													// change. Icarus, 50(2-3), pp.381-407.
												

	float tCO2 <- 0.0; //update: 0.004 * ((pCO2)* Pa2bar / ga )^0.4551; // Marinova et.al 2005 [Pa = kg*m-2 * N*kg-1 
															// CO2 equivalent grey opacity 						
	float tau <- 0.0; // update: tCO2 + tauDust;
	
	float Ts;									// calculated mars surface temperature
	float Tps <- 0.0; //update: Ts + Gamma_HT * height;	// potential temperature
	float Tsat <- 0.0;
	

	float prevTs <- 0.0;								// to keep the temperature from the previous step
	float eddy_Ts <- 0.0;
	float hadley_Ts <- 0.0;
	
	float ground_delta_En <- 0.0;
	float subground_Ts <- 0.0;
		
	float height_diff <- 0.0;
	float hadley_ht <- 0.0;
	
	float insol;
	float tauDustN;
	
	float kat;
	float zachSlonca;
	float OMEGA;
	
	float co2_excess <- 0.0;
	float latentCO2 <- 0.0;
	
	float insol_en;
	float rad_en;
	float green_en ;
	float prev_en;
	
	list<int> neigh <- list_with(6,0);
	list<int> azim <- list_with(6,0); 
			
	aspect base_co2{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb(255, (co2_column - 80)/2.0, 255); // border: #black;	
	}
	aspect base_Ts{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb((Ts - 80), 255, 255); // border: #black;	
	}
	aspect plain_co2{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(0, (co2_column - 80)/2.0, 0) border: #black;	
	}
	aspect plain_Ts{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: id_cell = selected_cell?#green:(frozenCO2>1e-8?#white:
			(Ts >= 273?#green:rgb((Ts - 120), 0, 0))
			) border: #black;	
	}
	aspect plain_energy{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(energy/3.0, 0, -energy/3.0) border: #black;	
	}
	aspect plain_FrozenCo2{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} 
			color: frozenCO2>1e-8?rgb(128+frozenCO2*10000, 128+frozenCO2*10000, 128+frozenCO2*10000):#black border: #black;	
	}
	aspect sphere_frozenCo2{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: frozenCO2>1e-8?rgb(128+frozenCO2*10000, 128+frozenCO2*10000, 128+frozenCO2*10000):#black; // border: #black;	
	}
	
	
	
	reflex initial when: cycle = 0 {
		neighbours <- cell where (each.id_cell in neigh);
		Ts <- temp;
		prevTs <- temp;
	}
	reflex finite_element when: cycle > 0 and cycle mod 2 = 1
	{
		float param <- (1-spec) * 4.0/6.0 + spec * 4.0/5.0; // jezeli spec to 4/5 inaczej 2/3
		float neigh_Ts_sum <- 0;
		loop n over: neighbours {
			neigh_Ts_sum <- neigh_Ts_sum + n.Tps;
		}
		eddy_Ts <- K_* (param * neigh_Ts_sum - 4 * Tps );
					 //-(beta_dw / delta_h)*(neighbours[beta_az].Tps - Tps)  
					//-(alfa_dw / delta_h)*(neighbours[alfa_az].Tps - Tps) 
					//+ K_* (param * neigh_Ts_sum - 4 * Tps );
	}
	reflex hadley_diff_cnt 
		when: cycle > 0 and cycle mod 2 = 1 
				and hadleyModel = 1 
		{
		hadley_Ts <- had_al[id_cell-1, sol_year];
	}
	
	/**
     * w modelu zaprooponowanym przez Chrisa wszystkie fluxy są per jednostka powierzchni!
     */
	reflex _balance when: cycle > 1 and cycle mod 2 = 0 {
		heat_flux <- energy;
		
		prevTs <- Ts;	
		
	   	float fi <- latitude;
		if (fi = 90.0)  { fi <- 89.99; } 
		if (fi = -90.0) { fi <- -89.99; }
		
		kat <- nachylenieOsi * sin(sol_lon);
		zachSlonca <- -1.0 * tan(fi) * tan(kat);
		OMEGA <- 0.0;
		
		float liczbaGodzSlonecznych <- 0.0;
		if (zachSlonca >= 1.0) {
			liczbaGodzSlonecznych <- 0.0;
			OMEGA <- 0.0;	
		} else if (zachSlonca <= -1) {
			liczbaGodzSlonecznych <- 24.0;
			OMEGA <- 180.0;
		} else {
			OMEGA <- acos(zachSlonca);
			liczbaGodzSlonecznych <- OMEGA * 2.0/15.0;
		}
		
		insol <- (1.0/#pi) * S0 * solar_const_modifier * ((1 + e * cos((sol_lon - 248) mod 360))^2) / ((1 - e^2)^2) 
					* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sol_lon) 
					* OMEGA * 2 * #pi / 360.0 + cos(fi) * cos(kat) * sin(OMEGA)
					) ; 		
		
		prevTs <- Ts;
		bool pCO2greater <- false;
		bool check <- true;
		float co2_excess_factor <- 1.0;
		tauDustN <- tauDust;
		
		loop times: 10 {
			
			pCO2 <- atmCO2 * exp(- height / scaleHeight ); // 0.006 * exp(- height / scaleHeight); // Pressure in [bar];
			if (pCO2 < 1e-8) {pCO2 <- 1e-8; }
			PsatCO2 <- 1.2264e7 * exp(-3167.8 / Ts); 
			Tsat <- -3167.8 / ln (pCO2 / 1.2264e7); // temperatura nasyceina - jaka powinna być wzgledem ciśnienia
			
			if (variableTauDust){
				tauDustN <-  tauDust * pCO2 / atmCO2;	
			}		
			
			//tCO2 <- 0.004 * ((pCO2)* Pa2bar / ga )^0.4551;
			tCO2 <- GHG_A * ((pCO2)* Pa2bar / ga )^GHG_B;
			tau <- tCO2 + tauDustN;
		
			latentCO2 <- 0.0;
			
			if (albedoFile){ // czytaj albedo z pliku
				if (sunGlass_Ice = 0) { // sunGlass_Ice == 0
					insol_en <- (1.0 - dustSunGlasses) *( (frozenCO2 > 1e-8)?(1 - (albedoIce - albedoGlobal + albedo)):(1 - albedo ) )* insol;  				// insolation ENERGY	
				} else { // sunGlass_Ice == 1
					insol_en <- (1.0 - ( (frozenCO2 > 1e-8)?dustSunGlasses:0.0)) *( (frozenCO2 > 1e-8)?(1 - (albedoIce - albedoGlobal + albedo)):(1 - albedo ) )* insol;  				// insolation ENERGY
				}
			} else { //przyjmij albedo stałe
				if (sunGlass_Ice = 0) { // sunGlass_Ice == 0
					insol_en <- (1.0 - dustSunGlasses) *( (frozenCO2 > 1e-8)?(1 - albedoIce):(1 - albedoGlobal ) )* insol;  				// insolation ENERGY	
				} else { // sunGlass_Ice == 1
					insol_en <- (1.0 - ( (frozenCO2 > 1e-8)?dustSunGlasses:0.0)) *( (frozenCO2 > 1e-8)?(1 - albedoIce):(1 - albedoGlobal ) )* insol;  				// insolation ENERGY
				}
					
			}
			insol_en <- insol_en * (1.0 - total_sg);
			//insol_en <- (1 - albedoGlobal ) * insol;  				// insolation ENERGY
			//insol_en <- (1 - albedo) * insol;  				// insolation ENERGY
			
			rad_en   <- emissivity * sigma             * prevTs^4;										// radiation of the planet
			//                       [W * m-2 * K-4 ]  * [K4]
			//						 [W * m-2]
			
			//green_en <- sigma * prevTs^4 * (1 - exp(-tau)) / (1 + 0.75 * tau);				// greenhouse effect
			green_en <- sigma * prevTs^4 * ( (7/8 + 3/2 * tau) * (1 - exp(-2*tau)) - 3/4 * tau  ) / (1 + 3/4*tau) ;	// greenhouse effect - mejl CHrisa 2023-03-13
			
			prev_en  <- prevTs * (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity) / delta_t; 
			// 1000 - 1500 Pa - exchangable CO2
			// total CO2 in atmosphere - atmoshpere + condension + polar caps
			// every cell has atmoshpere relating to pTotal
			// pres of hex <- pTotal 
			
			//if (id_cell = 1788 or id_cell = 333) {
			//	write "<"+id_cell+"> ["+n+"] pCO2 = " + pCO2 + ", pSat = " + PsatCO2 + ", Ts = " + Ts + ", Tsat = " + Tsat;
			//}
		
			
			
			if (pCO2 > PsatCO2){ // zbyt dużo CO2 - należy go "zabrać" z atmosfery i dodać do czapy lodowej
										// energia hexa maleje
				pCO2greater <- true;							
				
				//co2_excess <- totalCO2 * CO2freezingCoeff * ln(Tsat/Ts); // [bar] 
				//co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * (exp(exp_param * Tsat/Ts) - exp(exp_param))  ;
				if (latitude < sdsg_lat){
					co2_excess <- co2_excess_factor * totalCO2 * spCO2defreze * ( Tsat - Ts )  ;
				} else {
					co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * ( Tsat - Ts )  ;	
				}
				
				
				// Arhenius equation: A~10^10, Ea~30kJ/mol, MolMass = 44.0095
				//co2_excess <- CO2freezingCoeff * 10^10 * exp(-30000/(8.314462 * Ts))*44.0095/1000 / (marsArea / 4002);
				
				
				if (co2_excess > (atmCO2/4002)){
					co2_excess <- atmCO2/4002;
				}
				
				frozenCO2 <- frozenCO2 + co2_excess; // [bar]
				
				//resublimateCO2 <- co2_excess * 733932 * delta_t; // energia sublimacji [Wm-2] - bez mnozenia przez delta_h2 - energia per jednostka pow.
				// energia rośnie
				latentCO2 <- latentCO2 + co2_excess * 613000 * 4002;
//								  kg*m-2 * J/s * kg-1 = J/s*m-2
// dzielenie przez delta_t żeby dostać wyniki w Watach 
// The heat of sublimation for carbon dioxide - 32.3 kJ/mol.
// masa molowa CO2 - 44.0095 g/mol
// The heat of sublimation for carbon dioxide = 32.3 kJ/mol * (1/44.0095) mol/g = 0,73393 kJ/g = 733,93 kJ/kg = 733932 J/kg

				//if (id_cell = 1788 or id_cell = 333) {
				//	write "<"+id_cell+"> pCO2 > PsatCO2, co2_excess = " + co2_excess + ", latent en = " + latentCO2 + ", frozen = " + frozenCO2;
				//	
				//}
				
						
			} else if (pCO2  < PsatCO2){ // za mało CO2 - jeśli jest czapa lodowa - zabieramy trochę CO2 z czapy
			
				pCO2greater <- false;
				
				//co2_excess <- totalCO2 * CO2freezingCoeff * ln (Ts/Tsat); // [bar]
				//co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * (exp(exp_param * Ts/Tsat) - exp(exp_param)) ;
				//co2_excess <- CO2freezingCoeff * 10^10 * exp(-30000/(8.314462 * Ts))*44.0095/1000/ (marsArea / 4002);
				
				if (latitude < sdsg_lat){
					co2_excess <- co2_excess_factor * totalCO2 * spCO2defreze * ( Ts - Tsat ) ;
				} else {
					co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * ( Ts - Tsat ) ;
				}
				
				
				if (co2_excess > frozenCO2){
					co2_excess <- frozenCO2;
				}
				
				frozenCO2 <- frozenCO2 - co2_excess; // w [bar] !!!
				
				// energia maleje			
				latentCO2 <- latentCO2 - co2_excess * 613000 * 4002;				

				//if (id_cell = 1788 or id_cell = 333) {
				//	write "<"+id_cell+"> pCO2 < PsatCO2, co2_excess = " + (-co2_excess) + ", latent en = " + latentCO2 + ", frozen = " + frozenCO2;
				//	
				//}

			}
			
			if (latitude <= kreg_lat and latitude >= -kreg_lat){
				ground_delta_En <- k_reg   * (prevTs - subground_Ts) / deepSoil;
			} else {
				ground_delta_En <- k_reg * kreg_mod * (prevTs - subground_Ts) / deepSoil;
			}
			
		 	if (deepGroundModel = 0) { // jeśli tryb bez modelu głębokiego - zreruj
		 		ground_delta_En <- 0.0;
		 	}	//W m-2			// W m-1 K-1      W m-2                          m-2 = W2 m-5 K-1  
			    // kg-1 m3 W-1 s-1 kg K * W2 m-5 K-1 = W m-2
			
			energy <- insol_en - rad_en + green_en + eddy_Ts + hadley_Ts  + latentCO2 - ground_delta_En;  // latent heat
//@			energy <- insol_en - rad_en + green_en + eddy_Ts + hadley_Ts  + latentCO2;  // latent heat
			
									   
			//                        K  =  W m-2            s  m-1 kg -1 m3 W-1 s-1 kg K
			//					      K  =  K
				   
			float tmpTs <- prevTs			// previous step temp. 
				+  energy * delta_t / (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity); // mass m-2 * conduction - without multiplication by area - its flux
//										  Wm-2 / J*kg-1*K-1 *  kg * m-2          + J*kg-1*K-1 * m         * kg * m-3
//										  Wm-2 1/ J*m-2*K-1                       + J*K-1*m-2
//										  Wm-2 /K-1 * m2 * J-1
//										  Ks-1 
			if (tmpTs < constrTemp) {
				tmpTs <- constrTemp;
			}
		
			// check if co2_excess is too large
			float tmpPsatCO2 <- 1.2264e7 * exp(-3167.8 / tmpTs);
			float tmpPCO2 <- (atmCO2 + co2_excess * (pCO2greater?(1.0):(-1.0)) ) * exp(- height / scaleHeight);
			
			if ( (tmpPCO2 > tmpPsatCO2 and pCO2greater) or (tmpPCO2 < tmpPsatCO2 and !pCO2greater)) {break; }
			//if (pCO2greater or (tmpPCO2 < tmpPsatCO2 and !pCO2greater)) {break; }
			//if (abs(tmpPCO2 - tmpPsatCO2 ) < 1e-2) { break; }
			else {
				co2_excess_factor <- co2_excess_factor * 0.5; 
			}
		} 
		Ts <- prevTs			// previous step temp. 
				+  energy * delta_t / (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity); // mass m-2 * conduction - without multiplication by area - its flux
//										  Wm-2 / J*kg-1*K-1 *  kg * m-2          + J*kg-1*K-1 * m         * kg * m-3
//										  Wm-2 1/ J*m-2*K-1                       + J*K-1*m-2
//										  Wm-2 /K-1 * m2 * J-1
//										  Ks-1 
		if (Ts < constrTemp) {
			Ts <- constrTemp;
		}
				
		//subground_Ts <- subground_Ts + ground_delta_En * delta_t / (deepSoil * cSoil * soilDensity);		
		//subground_Ts <- subground_Ts + ground_delta_En / (11); 10a		
		subground_Ts <- subground_Ts + delta_t * ground_delta_En / (cSoil * soilDensity * deepSoil);
		// K                  K + s W m-2 / (J*kg-1*K-1 * m         * kg * m-3)
		// K                  K + s W m-2 / (K-1 m-2  W s )
		// K                  K + s W m-2 K m2  W-1 s-1 )
		// K                  K +  K  
		

		Tps <- Ts - Gamma_HT * height; // poprawka - zmiana znaku na "-" 2023-05-17
							
		heat_flux <- heat_flux - energy;	// obliczenie przepływu - różnicowo, pewnie trzeba dodać sol * delta_h, ale to później
		
			
	}

	reflex diff_v1 when: log_calib and id_cell = 333 and ( (sol_year >= 180 and sol_year <= 286) or (sol_year >= 303 and sol_year <= 553)){
		//Vik 1
		int offset <- 0;
		if      (sol_year <= 286) {offset <- 180;}
		else if (sol_year <= 553) {offset <- 197;}
		
		mape_tv1 <- mape_tv1 + abs(vik1[1, sol_year - offset] - Ts)/vik1[1, sol_year - offset];
		mape_pv1 <- mape_pv1 + abs(vik1[2, sol_year - offset]/1e5 - pCO2)/(vik1[2, sol_year - offset]/1e5);		
		
		rmse_tv1 <- rmse_tv1 + (vik1[1, sol_year - offset] - Ts)^2;
		rmse_pv1 <- rmse_pv1 + (vik1[2, sol_year - offset]/1e5 - pCO2)^2;
		
				
	}
	
	reflex diff_v2 when: log_calib and id_cell = 3090 and ( (sol_year >= 219 and sol_year <= 346) 
										   or (sol_year >= 354 and sol_year <= 357)
										   or (sol_year >= 362 and sol_year <= 596)
	){
		//Vik 2
		int offset <- 0;
		if      (sol_year <= 346) {offset <- 219;}
		else if (sol_year <= 357) {offset <- 219+8;}
		else if (sol_year <= 596) {offset <- 219+8+5;}
		
		mape_tv2 <- mape_tv2 + abs(vik2[1, sol_year - offset] - Ts)/vik2[1, sol_year - offset];
		mape_pv2 <- mape_pv2 + abs(vik2[2, sol_year - offset]/1e5 - pCO2)/vik2[2, sol_year - offset]/1e5;		
				
		rmse_tv2 <- rmse_tv2 + (vik2[1, sol_year - offset] - Ts)^2;
		rmse_pv2 <- rmse_pv2 + (vik2[2, sol_year - offset]/1e5 - pCO2)^2;
	}
	
	
	
	
	reflex f when: id_cell = selected_cell {
		tempFrozen <- frozenCO2;
		tempTemp <- Ts;
		tempTsat <- Tsat;
		tempCO2 <- pCO2; 
		tempSat <- PsatCO2;
		tempEnergy  <- energy;
		
		eddyEn <- eddy_Ts;
		hadEn <- hadley_Ts;
		latentEn <- latentCO2;
		greenEn <- green_en;
		radEn <- - rad_en;
		insolEn <- insol_en;
		prevEn <- prev_en;
		
		grEn <- ground_delta_En;
		sgrTs <- subground_Ts;	
	}
	
	reflex update_mtx {
		if (year = 0) {
			temp_mtx[0,sol_year] <- Ts;
			pres_mtx[0,sol_year] <- pCO2;
		} else if (year = 1) {
			temp_mtx[1, sol_year] <- Ts;
			pres_mtx[1, sol_year] <- pCO2;
		} else {
			temp_mtx[0, sol_year] <- temp_mtx[1, sol_year];
			pres_mtx[0, sol_year] <- pres_mtx[1, sol_year];
			temp_mtx[1, sol_year] <- Ts;
			pres_mtx[1, sol_year] <- pCO2;
		}
	}
	reflex check_diff when: log_calib and year > 0 and sol_year = 667 {
		mtx_diff <- 0.0;
		loop s from: 0 to: 667 {
			mtx_diff <- mtx_diff + (temp_mtx[1,s] - temp_mtx[0,s])^2/668 + (pres_mtx[1,s] - pres_mtx[0,s])^2/668;
		}
	}

}

experiment main_experiment until: (cycle > 6680) 
{
	
	parameter "Start sollon" var: sollon_number min: 0 max: 45;
	parameter "Log every sol" var: log_every_sol;
	parameter "Log from the year no" var: log_first_year;
	parameter "Selected cell" var: selected_cell;
	parameter "Selected cells" var: selected_cells;
	parameter "Sublimation paramter" var: CO2freezingCoeff;
	parameter "Sublimation paramter for south pole area" var: spCO2defreze;
	parameter "Eddy heat transfer coeef" var: Kcoeff;
	parameter "Hadley heat transfer coeef" var: hadleyParam;
	parameter "Obliquity" var: nachylenieOsi;
	parameter "Gamma (Temperature - height)" var: Gamma_HT;
	parameter "Const temp [K]" var: constrTemp;
	parameter "Soil thickness" var: soilDepth;
	parameter "Subsurface deep" var: deepSoil;
	parameter "Subsurface initial temp" var: initDeepGrdTemp;
//	parameter "Deep soil mass modificator" var: soilMassMod;
	parameter "Hudely model, 0-simple, 1-Alison" var: hadleyModel;
	parameter "Deep ground model, 0-none, 1-exists" var: deepGroundModel;
	parameter "K regolith modifier for areas near the poles" var: kreg_mod;
	parameter "Boundary lattitudes for K regolith modifier" var: kreg_lat;
	parameter "Scale height" var: scaleHeight;
	

	parameter "Soil emissivity" var: emissivity;
	parameter "Total CO2 on Mars" var: totalCO2;
	parameter "Soil gray opacity" var: tauDust;
	parameter "Albedo" var: albedoGlobal;
	parameter "Ice albedo" var: albedoIce;
	parameter "albedo from file" var: albedoFile;
	parameter "Total sunglasses" var: total_sg;
	parameter "Put sunglasses when ice only" var: sunGlass_Ice;
	parameter "Dust subglasses on south pole" var: south_dsg;
	parameter "Dust subglasses on north pole" var: north_dsg;
	parameter "North latitude of pole" var: ndsg_lat;
	parameter "South latitude of pole" var: sdsg_lat;
	parameter "Constant (true) or ice-dependent (false) dust subglasses" var: sun_glass_param;
	parameter "Variable Tau Dust: 0 - constant, 1 - variable" var: variableTauDust;
	parameter "solar constant modifier" var: solar_const_modifier;
	parameter "solar pole glass mode" var: sun_glass_mode;
	parameter "GHG A" var: GHG_A;
	parameter "GHG B" var: GHG_B;
	
	parameter "Folder for result" var: outdir;	
	parameter "Log output" var: log_output;
	parameter "Log to CSV" var: log_csv;
	parameter "Log Viking cells" var: log_vik;
	parameter "Log calibration" var: log_calib;
	 
	
	output {
	
		display mars_plain_Ts type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_Ts;
		}	
		
	
		display chart2 refresh: every(2#cycles) {
			chart "Temperature" type: series background: #white style: exploded {
				data "Surface temperature" value: tempTemp color: #red;  
				data "Saturation temperature" value: tempTsat color: #green;  
				data "Subground Ts" value: sgrTs color: #black;
				data "Flux" value: grEn color: #blue;
			}
		}
		display chart3 refresh: every(2#cycles) {
			chart "CO2 [bar]" type: series background: #white style: exploded {
				data "CO2" value: tempCO2 color: #black;  
				//data "Sat press" value: tempSat color: #red;
				data "Frozen CO2" value: tempFrozen color: #blue;  
			}
		}
		display chart4 refresh: every(2#cycles) {
			chart "Energy" type: series background: #white style: exploded {
				data "Total Energy" value: tempEnergy color: #black;  
				data "Latent Energy" value: latentEn color: #blue;
				data "Greenhouse Energy" value: greenEn color: #green;
				data "Radiation Energy" value: radEn color: #red;
				data "Insolation Energy" value: insolEn color: #orange;
//				data "Preserved Energy" value: prevEn color: #gray;
				data "Eddy Energy"  value: eddyEn color: #purple;
				data "Hadley Energy" value: hadEn color: #cyan;
				data "Ground Energy" value: grEn color: #orange;			
			}
		}
		/*display chart6 refresh: every(2#cycles) {
			chart "Total CO2" type: series background: #white style: exploded {
				data "Total CO2 mass" value: totalCO2 * delta_h2 * Pa2bar / ga  color: #red;  
				data "Total frozen CO2 mass" value: sumFrozen * delta_h2 * Pa2bar / ga color: #blue;  
				data "Atmosphere CO2 mass" value: atmCO2 * delta_h2 * Pa2bar / ga color: #yellow;  
			}
		}*/
		display chart6 refresh: every(2#cycles) {
			chart "Total CO2 [bar]" type: series background: #white style: exploded {
		//		data "Total CO2" value: totalCO2   color: #red;  
				data "Total frozen CO2" value: sumFrozen color: #blue;  
				data "Atmosphere CO2" value: cell mean_of(each.pCO2) color: #yellow;
				data "MCD CO2" value: calib[4, sol_year]/1e5 color: #gray;  
			}
		}
		display chart6a refresh: every(2#cycles) {
			chart "frozen / CO2" type: series background: #white style: exploded {
				data "Percent of frozen CO2" value: sumFrozen / totalCO2 color: #blue;  
			}
		}
		
		
		/*display chart7 refresh: every(2#cycles) {
			chart "Total Energy [10e15 J]" type: series background: #white style: exploded {
				data "Total Energy [10e15 J]" value: cell sum_of(each.energy) color: #red;  
			}
		}*/
		display chart8 refresh: every(2#cycles) {
			chart "Mean temp" type: series background: #white style: exploded {
				data "Mean temp" value: cell mean_of(each.Ts) color: #red;  
				data "MCD temp" value: calib[1, sol_year] color: #gray;
				data "Minimal temp" value: cell min_of(each.Ts) color: #orange;
 			}
		}
		/*display chart9 refresh: every(2#cycles) {
			chart "Mean insol" type: series background: #white style: exploded {
				data "Mean insol" value: cell mean_of(each.insol) color: #red;  
			}
		}*/
		
		monitor "Mean MARS greenhouse temperature"                   name: mean_greenhouse_temp value: cell mean_of(each.Ts);
		monitor "Varianve of MARS greenhouse temperature"            name: var_greenhouse_temp value: cell variance_of(each.Ts);
		monitor "Max MARS greenhouse temperature"            		 name: max_greenhouse_temp value: cell max_of(each.Ts);
		monitor "Min MARS greenhouse temperature"            		 name: min_greenhouse_temp value: cell min_of(each.Ts);
		
		monitor "Mean MARS potential temperature"                    name: mean_potential_temp value: cell mean_of(each.Tps);
		
		monitor "Mean MARS CO2 pressure"                   			 name: mean_pCO2 value: cell mean_of(each.pCO2);
		monitor "Varianve of MARS CO2 pressure"            			 name: var_pCO2 value: cell variance_of(each.pCO2);
		monitor "Max MARS CO2 pressure"            		 			 name: max_pCO2 value: cell max_of(each.pCO2);
		monitor "Min MARS CO2 pressure"            		 			 name: min_pCO2 value: cell min_of(each.pCO2);
		
		monitor "Sum of frozen CO2" 								 name: sum_frozenCO2 value: cell sum_of(each.frozenCO2 * delta_h2) * Pa2bar / ga;
		monitor "Sum of energy" 								 	 name: sum_energy value: cell sum_of(each.energy);
		
		
		monitor "Hexes with frozen CO2"            		 			 name: no_hex_frozenCO2 value: cell count(each.frozenCO2 > 1e-8);
		monitor "Number of hexes biologically habitable (T > -25 C)" name: biol_habitable_hex value: cell count(each.Ts > 248.15);
		monitor "Number of hexes with unfreezed water (T > 0 C)"     name: biosphere_hex      value: cell count(each.Ts > 273.15);
		
		monitor "Mean equatorial temperature"						 name: mean_equatorial_temp value: equ_list mean_of (each.Ts);
		monitor "Mean north pole temperature"						 name: mean_Npole_temp value: N_poles_list  mean_of (each.Ts);
		monitor "Mean south pole temperature"						 name: mean_Spole_temp value: S_poles_list  mean_of (each.Ts);
		monitor "Mean +45 temperature"						 		 name: mean_45p_temp value: _45p_poles_list  mean_of (each.Ts);
		monitor "Mean -45 temperature"						 		 name: mean_45m_temp value: _45m_poles_list  mean_of (each.Ts);
		
		monitor "Vik 1 temperature"									 name: vik1_temp value: cell[333].Ts; 
		monitor "Vik 2 temperature"									 name: vik2_temp value: cell[3090].Ts; 
		monitor "Vik 1 pressure"									 name: vik1_pres value: cell[333].pCO2; 
		monitor "Vik 2 pressure"									 name: vik2_pres value: cell[3090].pCO2; 
		monitor "Vik 1 gnd flux"									 name: vik1_gflux value: cell[333].ground_delta_En; 
		monitor "Vik 2 gnd flux"									 name: vik2_gflux value: cell[3090].ground_delta_En; 
		
		monitor "Daily insolation"									 name: daily_insol value: cell mean_of(each.insol);
//		monitor "Fitness"											 name: fit value: fitness;
//		monitor "Year to year difference (RMSE)"					 name: yty value: yty_diff;
		
				
		monitor "Sols"     name: solno value: cycle / 2 / step_sol ;
//		monitor "HN"     name: hn value: hadleyN;
//		monitor "HS"     name: hs value: hadleyS;
		
		monitor "Temp in selected cell"     name: tempSelCell value: tempTemp;
		monitor "CO2 Press in selected cell"     name: presSelCell value: tempCO2;
		monitor "CO2 frozen in selected cell"     name: frozSelCell value: tempFrozen;
		monitor "Latent heat in selected cell"     name: latentSelCell value: latentEn;
		monitor "Insolation heat in selected cell"     name: insolSelCell value: insolEn;
		monitor "Ground delta Energy"     name: grEn value: cell mean_of (each.ground_delta_En);
		monitor "Subground Energy"     name: sgrTs value: cell mean_of (each.subground_Ts);
		
	}
}
/*
experiment hillclimb type: batch keep_seed: true repeat: 1 until: ( cycle >= 6680 ) {
    parameter 'Hadley heat transfer coeef' var: hadleyParam min: 0.0000000005 max: 0.0000000095 step: 0.0000000005;
    parameter "co2 total" var: totalCO2 min: 0.006 max: 0.008 step: 0.0001;
    parameter "freezing coeff" var: CO2freezingCoeff min: 0.0000005 max: 0.0000095 step: 0.0000005;
    parameter "log calib" var: log_calib;
    
	
    method hill_climbing iter_max: 50  minimize: fitness  ; 
    
    output{
    	monitor "Fitness"     name: fit value: fitness;
    	
    	display chart  {
			chart "Fitness" type: series background: #white style: exploded {
				data "Fitness" value: fitness color: #red;  
			}
		}
    }
}  */
