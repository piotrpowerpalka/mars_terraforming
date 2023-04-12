/**
* Name: MarsEnergyModel
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
* 

*/


model MarsEnergyModelBucket


/* Insert your model definition here */

global {
	float f <- 10e-4; 			// [sec-1] - Coriolis parameter
	float R <- 0.19e-3; 		// [J kg-1 K-1] gas constant = 0.19×107 ergs gm-1 K-1
	float S <- 2e-7;			// [K *m-1] static stability
	float L <- 5.3e6;			// [m] equator-pole distance
	float ro <- 2.4733;		// [kg m-3] CO2 density in 220 K
	
	float totalCO2 <- 0.0067; // [bar] total amount of CO2 on Mars
	float sumFrozen update: cell sum_of(each.frozenCO2); 
	float atmCO2 update: totalCO2 - sumFrozen;  // [bar] mean CO2 pressure
	
	int aspect_mode <- 1;				// tryb wyświetlania
	int sollon_number <- 0;
    int model_number <- 3;
    int selected_cell <- 333;
    list<int> selected_cells <- [1788, 2859];
    
    int numOfHexes <- 4002;
    float Tincrease <- 0.0;
    
    int GHGregress <- 0;
    float Gamma_HT <- 10/1000; // Kelvin / m - potential temperature parameter
    //float Gamma_HT <- ga/cCO2; // 0.004751 wzór wg. Handbook on Atmospheric Diffusion
     
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
	
	float tauDust <- 0.0;
	float albedoGlobal <- 0.2;	
		
	string outdir <- "../results/";
	    
	int log_every_sol <- 1;
	bool log_output <- false;
	bool log_csv <- false;
	bool log_vik <- false;
	
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
	float soilDepth <- 0.000075; // [m]
	float soilDensity <- 1600; // [kg m-3]
	float emissivity <- 0.9; // [W m-2]
	
	//float insolParam  <- 1.0; // parameter limiting insoloation - to decrease temp
	float insolParam  <- 0.55; // parameter limiting insoloation - to decrease temp
	
	float delta_t <- 88775.0;		// 1 martian sol [s]
	float delta_h2 <-  marsArea / numOfHexes;	// area of single hex
	float delta_h <- delta_h2 / #pi;			// radius of single hex (approx.)
	float Kcoeff <- 10000.0;
	float K_ <- Kcoeff * delta_t / delta_h2 / step_sol; // K * delta_t / delta_h2
	float Wind <- 0.0;
	
		
   	float hadleyN <- 0.0;
   	float hadleyS <- 0.0;
   	//float hadleyParam <- 5e-9;
   	float hadleyParam <- 5e-6;
   	
		
	int martianYear <- 668; // martian year [sol]
	float Pa2bar <- 1e5;			// pascal to bar
	
	// albedo
	file csv_albedo <- csv_file("../includes/heksy_albedo.csv",";",float,true);
	matrix<float> mat_albedo <- matrix<float>(csv_albedo.contents);
	
	file csv_hadley <- csv_file("../includes/hadley_heat_transport.csv",",",float,true);
	matrix<float> mat_hadley <- matrix<float>(csv_hadley.contents);
	
	
	int step_sol <- 1;	// number of iterations per sol
	int sol_year <- 0 update: (cycle/(2 * step_sol)) mod martianYear; // number of sol in martian year
    int sol_lon  <- 0 update: int(sol_year / 668.0 * 360.0); // recalculation sol no -> solar logitude
    
    
	float mean_greenhouse_temp 	<- 0.0 update: cell mean_of (each.Ts);  // mean greenhouse temperature
	float var_greenhouse_temp 	<- 0.0;
	float mean_equatorial_temp 	<- 0.0 update: equ_list mean_of (each.Ts);
	float mean_pole_temp 		<- 0.0 update: poles_list mean_of(each.Ts); 
	
	int biol_habitable_hex <- 0;		// number of hexes biologically habitable (T > -25 C)
	int biosphere_hex <- 0; 			// number of hexes unfreezed water (T > 0 C)
	
    
	float CO2freezingCoeff <- 0.0000025 / step_sol; // percentage part of CO2 being sublimated / resublimated per sol
    float CO2defreezingCoeff <- 0.0000025 / step_sol; // percentage part of CO2 being resublimated per sol
    
    list<cell> poles_list <- nil;
    list<cell> equ_list <- nil;
 
	
    init {
			
		
	
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
				hadley_ht <- mat_hadley[id_cell-1];
								
  				location <- {longitude, latitude};
    			neigh <- [n1, n2, n3, n4, n5, n6];
				azim  <- [a1, a2, a3, a4, a5, a6];
				Ts <- temp;
				
				if (spec = 1) {
					remove index: 5 from: neigh;
					remove index: 5 from: azim;
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
    		
    	poles_list <- cell where (each.latitude >= 90 or each.latitude <= -90);
		equ_list   <- cell where (each.latitude >= -2 and each.latitude <= 2);

    }
    
    /**
     * save data every year
     */
    reflex saving_output when: log_output and (
    	cycle mod (2 * step_sol * log_every_sol) = 1 // or cycle mod 24 = 1
    	//cycle mod 56 = 1
    	) {    	
    	loop n over: cell   {
    		int sol <-  (cycle -1) / (2 * step_sol);
 	 			write "" + ( (cycle -1) / (2 * step_sol) ) + ";" + n.id_cell+ ";" + 
	 						n.Ts + ";" + (n.pCO2 - sumFrozen/numOfHexes) + ";"  +
    	          	 		n.frozenCO2 + ";" + n.heat_flux + ";" + 
    	          	 		(sol mod 668) + ";" + (sol - (sol mod 668))/668 + ";" +
    	          	 		n.insol
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
	 						n.Ts + ";" + (n.pCO2 - sumFrozen/numOfHexes) + ";"  +
    	          	 		n.frozenCO2 + ";" + n.heat_flux  + ";" + 
    	          	 		(sol mod 668) + ";" + (sol - (sol mod 668))/668 + ";" +
    	          	 		n.insol
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
	    	    	n.Ts, n.insol, n.energy, (n.pCO2 - sumFrozen/numOfHexes), n.frozenCO2, n.heat_flux
	    	    ] to: outdir + "/year_" + int( cycle / 1336 )  + "_sol_" + (cycle / 2) mod 668 + "_cell_" + n.id_cell + ".csv" rewrite: true type: "csv";
    	    }    	   
     	}
    }
    
    
    reflex hadley_heat_transport when: cycle > 0 and cycle mod 2 = 1 {
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
    		ask n {n.hadley_Ts <- n.hadley_Ts * mean_greenhouse_temp^0.5 * (mean_equatorial_temp - mean_pole_temp)^2; }
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
		
	
	float longitude;
	float latitude;
	float height;
	float zonal_wind;
	float merid_wind;
	
	float atm_press;
	int spec;
	
	float co2_column;
	float next_step_co2_column;
 	float next_step_Ts;
 	float neigh_Ts_sum;
 	
 	float div_co2 <- 0.0;
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
	float pCO2 update: atmCO2 * exp(- height / 10000); // 0.006 * exp(- height / 10000); // Pressure in [bar]
	
	float PsatCO2 update: 1.2264e7 * exp(-3167.8 / Ts); //[bar] za Reference this from Eq 19 in Fanale et al. (1982) Fanale, F.P., Salvail, J.R., Banerdt, W.B. and 
													// Saunders, R.S., 1982. Mars: The regolith-atmosphere-cap system and climate 
													// change. Icarus, 50(2-3), pp.381-407.
												

	float tCO2 update: 0.004 * ((pCO2)* Pa2bar / ga )^0.4551; // Marinova et.al 2005 [Pa = kg*m-2 * N*kg-1 
															// CO2 equivalent grey opacity 						
	float tau update: tCO2 + tauDust;
	
	float Ts;									// calculated mars surface temperature
	float Tps update: Ts + Gamma_HT * height;	// potential temperature
	
	

	float prevTs <- 0.0;								// to keep the temperature from the previous step
	float eddy_Ts <- 0.0;
	float hadley_Ts <- 0.0;
	
	float height_diff <- 0.0;
	float hadley_ht <- 0.0;
	
	float insol;
		
	float kat;
	float zachSlonca;
	float OMEGA;
	
	float co2_excess <- 0.0;
	float resublimateCO2 <- 0.0;
	float sublimateCO2 <- 0.0;
	
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
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: id_cell = selected_cell?#green:(frozenCO2>0.001?#white:rgb((Ts - 80), 0, 0)) border: #black;	
	}
	aspect plain_energy{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(energy/3.0, 0, -energy/3.0) border: #black;	
	}
	aspect plain_FrozenCo2{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} 
			color: frozenCO2>0.001?rgb(128+frozenCO2*10000, 128+frozenCO2*10000, 128+frozenCO2*10000):#black border: #black;	
	}
	aspect sphere_frozenCo2{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: frozenCO2>0.001?rgb(128+frozenCO2*10000, 128+frozenCO2*10000, 128+frozenCO2*10000):#black; // border: #black;	
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
		
		insol <- insolParam * (24/#pi) * S0 * ((1 + e * cos((sol_lon - 248) mod 360))^2) / ((1 - e^2)^2) 
					* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sol_lon) 
					* OMEGA * 2 * #pi / 360.0 + cos(fi) * cos(kat) * sin(OMEGA)
					) / (24 ) ; 
		// [W * m-2]
	
		div_co2 <- 0.0;
		
		//insol_en <- (1 - ( (frozenCO2 > 0.001)?0.8:albedoGlobal) ) * insol;  				// insolation ENERGY
		insol_en <- (1 - albedoGlobal ) * insol;  				// insolation ENERGY
		//insol_en <- (1 - albedo) * insol;  				// insolation ENERGY
		
		rad_en   <- emissivity * sigma             * prevTs^4;										// radiation of the planet
		//                       [W * m-2 * K-4 ]  * [K4]
		//						 [W * m-2]
		
		//green_en <- sigma * prevTs^4 * (1 - exp(-tau)) / (1 + 0.75 * tau);				// greenhouse effect
		green_en <- sigma * prevTs^4 * ( (7/8 + 3/2 * tau) * (1 - exp(-2*tau)) - 3/4 * tau  ) / (1 + 3/4*tau) ;	// greenhouse effect - mejl CHrisa 2023-03-13
		
		
		prev_en  <- prevTs * (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity); 
		// 1000 - 1500 Pa - exchangable CO2
		// total CO2 in atmosphere - atmoshpere + condension + polar caps
		// every cell has atmoshpere relating to pTotal
		// pres of hex <- pTotal 
		
		energy <- insol_en - rad_en + green_en + eddy_Ts + hadley_Ts 
					- sublimateCO2 + resublimateCO2;  // latent heat
							   
		
		Ts <- prevTs			// previous step temp. 
				+  energy / (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity); // mass m-2 * conduction - without multiplication by area - its flux
//										  Wm-2 / J*kg-1*K-1 *  kg * m-2          + J*kg-1*K-1 * m         * kg * m-3
//										  Wm-2 1/ J*m-2*K-1                       + J*K-1*m-2
//										  Wm-2 /K-1 * m2 * J-1
//										  Ks-1 
							
		heat_flux <- heat_flux - energy;	// obliczenie przepływu - różnicowo, pewnie trzeba dodać sol * delta_h, ale to później
		
			sublimateCO2 <- 0.0;
			resublimateCO2 <- 0.0;
			
			if (pCO2 > PsatCO2){ // zbyt dużo CO2 - należy go "zabrać" z atmosfery i dodać do czapy lodowej
										// energia hexa maleje
				co2_excess <- totalCO2 * CO2freezingCoeff; // [bar] 
				
				if (co2_excess > atmCO2){
					co2_excess <- atmCO2;
				}
				
				div_co2 <- -co2_excess; // [bar]
				frozenCO2 <- frozenCO2 + co2_excess; // [bar]
				//resublimateCO2 <- co2_excess * delta_h2 * 289750; // energia sublimacji [J]
				resublimateCO2 <- co2_excess * 733932 / delta_t ; // energia sublimacji [Wm-2] - bez mnozenia przez delta_h2 - energia per jednostka pow.
//								  kg*m-2 * J/s * kg-1 = J/s*m-2
// dzielenie przez delta_t żeby dostać wyniki w Watach 
// The heat of sublimation for carbon dioxide - 32.3 kJ/mol.
// masa molowa CO2 - 44.0095 g/mol
// The heat of sublimation for carbon dioxide = 32.3 kJ/mol * (1/44.0095) mol/g = 0,73393 kJ/g = 733,93 kJ/kg = 733932 J/kg
						
			} else if (pCO2  < PsatCO2){ // za mało CO2 - jeśli jest czapa lodowa - zabieramy trochę CO2 z czapy
				co2_excess <- totalCO2 * CO2defreezingCoeff; // [bar]
				
				if (co2_excess > frozenCO2){
					co2_excess <- frozenCO2;
				}
				
				div_co2 <- co2_excess;
				frozenCO2 <- frozenCO2 - co2_excess; // w [bar] !!!
				//sublimateCO2 <- co2_excess * delta_h2 * 289750 ; // energia resublimacji [J]
				sublimateCO2 <- co2_excess * 733932 / delta_t; // energia resublimacji [Wm-2]	- bez mnozenia przez delta_h2 - energia per jednostka pow.			
								
			}
		
	}
	reflex f when: id_cell = selected_cell {
		tempFrozen <- frozenCO2;
		tempTemp <- Ts;
		tempCO2 <- pCO2; 
		tempSat <- PsatCO2;
		tempEnergy  <- energy;
		
		eddyEn <- eddy_Ts;
		hadEn <- hadley_Ts;
		latentEn <- - sublimateCO2 + resublimateCO2;
		greenEn <- green_en;
		radEn <- - rad_en;
		insolEn <- insol_en;
		prevEn <- prev_en;
	}
}

experiment main_experiment until: (cycle > 6680) 
{
	
	parameter "Start sollon" var: sollon_number min: 0 max: 45;
	parameter "Log every sol" var: log_every_sol;
	parameter "Selected cell" var: selected_cell;
	parameter "Selected cells" var: selected_cells;
	parameter "Sublimation paramter" var: CO2freezingCoeff;
	parameter "Resublimation paramter" var: CO2defreezingCoeff;
	parameter "Eddy heat transfer coeef" var: Kcoeff;
	parameter "Hadley heat transfer coeef" var: hadleyParam;
	parameter "Obliquity" var: nachylenieOsi;
	parameter "Gamma (Temperature - height)" var: Gamma_HT;
	parameter "Soil thickness" var: soilDepth;
	parameter "Soil emissivity" var: emissivity;
	parameter "Total CO2 on Mars" var: totalCO2;
	parameter "Soil gray opacity" var: tauDust;
	parameter "Insolation limitation" var: insolParam;
	parameter "Albedo" var: albedoGlobal;
	parameter "Folder for result" var: outdir;	
	parameter "Log output" var: log_output;
	parameter "Log to CSV" var: log_csv;
	parameter "Log Viking cells" var: log_vik;
	 
	
	output {
	
		display mars_plain_Ts type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_Ts;
		}	
		
	
		display chart2 refresh: every(2#cycles) {
			chart "Temperature" type: series background: #white style: exploded {
				data "Temperature" value: tempTemp color: #red;  
			}
		}
		display chart3 refresh: every(2#cycles) {
			chart "CO2" type: series background: #white style: exploded {
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
				//data "Preserved Energy" value: prevEn color: #gray;
				data "Eddy Energy"  value: eddyEn color: #purple;
				data "Hadley Energy" value: hadEn color: #cyan;
	
				
			}
		}
		display chart6 refresh: every(2#cycles) {
			chart "Total CO2" type: series background: #white style: exploded {
				data "Total CO2 mass" value: totalCO2 * delta_h2 * Pa2bar / ga  color: #red;  
				data "Total frozen CO2 mass" value: sumFrozen * delta_h2 * Pa2bar / ga color: #blue;  
				data "Atmosphere CO2 mass" value: atmCO2 * delta_h2 * Pa2bar / ga color: #yellow;  
			}
		}
		display chart6a refresh: every(2#cycles) {
			chart "frozen / CO2" type: series background: #white style: exploded {
				data "Percent of frozen CO2" value: sumFrozen / totalCO2 color: #blue;  
			}
		}
		
		
		display chart7 refresh: every(2#cycles) {
			chart "Total Energy [10e15 J]" type: series background: #white style: exploded {
				data "Total Energy [10e15 J]" value: cell sum_of(each.energy) color: #red;  
			}
		}
		display chart8 refresh: every(2#cycles) {
			chart "Mean temp" type: series background: #white style: exploded {
				data "Mean temp" value: cell mean_of(each.Ts) color: #red;  
			}
		}
			
		monitor "Mean MARS greenhouse temperature"                   name: mean_greenhouse_temp value: cell mean_of(each.Ts);
		monitor "Varianve of MARS greenhouse temperature"            name: var_greenhouse_temp value: cell variance_of(each.Ts);
		monitor "Max MARS greenhouse temperature"            		 name: max_greenhouse_temp value: cell max_of(each.Ts);
		monitor "Min MARS greenhouse temperature"            		 name: min_greenhouse_temp value: cell min_of(each.Ts);
		
		monitor "Mean MARS CO2 pressure"                   			 name: mean_pCO2 value: cell mean_of(each.pCO2 - each.frozenCO2);
		monitor "Varianve of MARS CO2 pressure"            			 name: var_pCO2 value: cell variance_of(each.pCO2 - each.frozenCO2);
		monitor "Max MARS CO2 pressure"            		 			 name: max_pCO2 value: cell max_of(each.pCO2 - each.frozenCO2);
		monitor "Min MARS CO2 pressure"            		 			 name: min_pCO2 value: cell min_of(each.pCO2 - each.frozenCO2);
		
		monitor "Sum of frozen CO2" 								 name: sum_frozenCO2 value: cell sum_of(each.frozenCO2 * delta_h2) * Pa2bar / ga;
		monitor "Sum of energy" 								 	 name: sum_energy value: cell sum_of(each.energy);
		
		
		monitor "Hexes with frozen CO2"            		 			 name: no_hex_frozenCO2 value: cell count(each.frozenCO2 > 1e-8);
		monitor "Number of hexes biologically habitable (T > -25 C)" name: biol_habitable_hex value: cell count(each.Ts > 248.15);
		monitor "Number of hexes with unfreezed water (T > 0 C)"     name: biosphere_hex      value: cell count(each.Ts > 273.15);
		
		monitor "Mean equatorial temperature"						 name: mean_equatorial_temp value: equ_list mean_of (each.Ts);
		monitor "Mean pole temperature"						 		 name: mean_pole_temp value: poles_list mean_of (each.Ts);
				
		monitor "Sols"     name: solno value: cycle / 2 / step_sol ;
		monitor "HN"     name: hn value: hadleyN;
		monitor "HS"     name: hs value: hadleyS;
		
		monitor "Temp in selected cell"     name: tempSelCell value: tempTemp;
		monitor "CO2 Press in selected cell"     name: presSelCell value: tempCO2;
		monitor "CO2 frozen in selected cell"     name: frozSelCell value: tempFrozen;
		
	}
}
