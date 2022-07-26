/**
* Name: MARSCAHregress
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
*/


model MARSChris

/* Insert your model definition here */

global {
	int aspect_mode <- 1;				// tryb wyświetlania
	int green_mode <- 5;
	int impact_model <- 1;
	int regolith_model <- 2;
    float moving_percentage_of_gas <- 0.01;
    float opening_angle <- 60.0;
	int sollon_number <- 0;
    int model_number <- 3;
    int height_diff_include <- 1;
    int green_diff_include <- 0;
    int regLimit <- 0;
    float regValue <- 1000.0; // [kgm-2] - value of CO2 reservoir for each hex
    
    int numOfHexes <- 4002;
    int logRegress <- 0;
    float Td <- 30.0; // Temperature increment to outgas 1/e of regolith
    float Tincrease <- 0.0;
    
    int GHGregress <- 0;
    float Ts_ene_coeff <- 100.0;
    float Gamma_HT <- 0.01; // Kelvin / m - potential temperature parameter 
	
	
    
    float prevTotalNH3 <- 0.0;
    float co2vik1 <- 0.0;
    float co2vik2 <- 0.0;
    float tempVik1 <- 0.0;
    float tempVik2 <- 0.0;
    float eneVik1 <- 0.0;
    float eneVik2 <- 0.0;    
    
    float eneCell <- 0.0;
    float co2Cell <- 0.0;
    float tempCell <- 0.0;
    
    string cells_affected <- ""; 				// cell affected by the effect
    list<int> cells_aff;
    list<int> effect_time;					// cycle of the efect happened
    list<int> effect_stop;
    
    float nh3_const_increase <- 0.0;		// pressure of nh3 increased every iter
    float ch4_const_increase <- 0.0;		// pressure of nh3 increased every iter
    float cfc_const_increase <- 0.0;		// pressure of nh3 increased every iter
    
    float nh3_abrupt_increase <- 0.0;		// pressure of nh3 increased every iter
    float ch4_abrupt_increase <- 0.0;		// pressure of nh3 increased every iter
    float cfc_abrupt_increase <- 0.0;		// pressure of nh3 increased every iter
    
	string outdir <- "../results/";
	
	float Rh <- 0.7; // Relative humidity
    float Rgas <- 8.314; // Gas constant for water
    float Lheat <- 43655.0; // Latent heat (?)
    float P0 <- 1.4E6; //Reference pressure (ext(19)?)
    
	float sumNH3;
	float sumCH4;
	float sumCFC;
	
	int log_every_sol <- 668;
	
	float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
    float sol_const <- 589.0; 		//Present day martian solar constant [Wm-2]
	float ga <- 3.72076; 			// Mars gravitional acceleration [ms-2] [Nkg-1]
	float AU <- 149597870700.0;		// astronomical unit [m]
	float sunTemperature <- 5780.0; 	// sun temperature [K]
	float sunRadius <- 695700000.0;	// sun radius [m]
	float marsRadius <- 3396200.0;	// mars radius [m]
	float marsArea <- 1.448e12;		// mars area [m2]\
	
	float eddyDiffParam <- 10; // Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_cfc <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_ch4 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_nh3 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_Ts <- 1000.0; // Eddy diffusion parameter [m2s-1] (approx.)  - temperature
	float eddy_co2_coeff <- 1000;
	
	
	float A_eddy <- 1.859e-3; //  // emprical coeff: https://en.wikipedia.org/wiki/Mass_diffusivity
	int eddy_calculated <- 0;
	
	float gl_k_co2 <- 0.0;
	float gl_k_nh3 <- 0.0;
	float gl_k_ch4 <- 0.0;
	float gl_k_cfc <- 0.0;
	
	float delta_t <- 88775.0;		// 1 martian sol [s]
	
	float delta_h2 <-  marsArea / numOfHexes;	// area of single hex
	float delta_h <- delta_h2 / #pi;			// radius of single hex (approx.)
	
	float K_ <- eddyDiffParam * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_CFC_ <- eddyDiffParam_cfc * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_CH4_ <- eddyDiffParam_ch4 * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_NH3_ <- eddyDiffParam_nh3 * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_Co2_Temp <- eddy_co2_coeff ; // * delta_t / delta_h2;
	
	float K_Ts <- eddyDiffParam_Ts * delta_t / delta_h2;
	
	
	int martianYear <- 668; // martian year [sol]
	float Pa2bar <- 1e5;			// pascal to bar
	
	// albedo
	file csv_albedo <- csv_file("../includes/heksy_albedo.csv",";",float,true);
	matrix<float> mat_albedo <- matrix<float>(csv_albedo.contents);
	
	// height regression parameters
    file csv_height_regression <- csv_file("../includes/height_regress.csv", ";", float, true);
    matrix<float> mat_hr <- matrix<float>(csv_height_regression.contents);
    
    file csv_solmars_dist <- csv_file("../includes/mars_sol_distance.csv", ";", float, true);
	matrix<float> mat_solmars_dist <- matrix<float>(csv_solmars_dist.contents);
	
	file csv_co2_mat <- csv_file("../includes/matrix-0-360-extvar_66.csv", ";", float, true);
	matrix<float> mat_co2 <- matrix<float>(csv_co2_mat.contents);
	    
    int sol_year <- 0 update: (cycle/2) mod martianYear;
    int sol_lon  <- 0 update: int(sol_year / 668.0 * 360.0);
    
    list<regression> co2_fct_mat <- list_with(668, nil);
	matrix<float> instances <- 0.0 as_matrix {2, numOfHexes}; 
	
	list<regression> co2_mult_fct_mat <- list_with(668, nil);
	matrix<float> mult_instances <- 0.0 as_matrix {2, numOfHexes}; 
	
	
	float mean_greenhouse_temp <- 0.0;  // mean greenhouse temperature
	float var_greenhouse_temp <- 0.0;
	int biol_habitable_hex <- 0;		// number of hexes biologically habitable (T > -25 C)
	int biosphere_hex <- 0; 			// number of hexes unfreezed water (T > 0 C)
	
	float regolith <- 1000  * ga / Pa2bar; // regolith CO2 inventory
    float C <- regolith * 0.006^(-0.275) * exp(149 / Td); //Normalization constant for regolith calculations, from McKay and Fogg
    
    float regolith_thickness <- 0.01; // [m] grubość warstwy regolitu, która nagrzewa się od słonca (0.001m = 1mm)
    float regolith_density <- 1600.0; // [kg/m^3] --  1.6 [g/cm^3] = 1.6 * 0.001 [kg] / 0.000001 [m^3] = 1600 [kg/m^3] za  
    //Behavior of Carbon Dioxide and Other Volatiles on Mars,  Robert B. Leighton and Bruce C. M
    float regHeatCap <- 3300; // regolith heat cappacity [J/kg * K], za
    //Behavior of Carbon Dioxide and Other Volatiles on Mars,  Robert B. Leighton and Bruce C. M
    
    float marsInclination <- 25.19; // [deg] 
    
	init {
		/**
		 * read arguments: 
		 * formats:
		 * 1528 -- a cell_id affected
		 * 1528;1224;1000 - cells affected
		 * 1528:668-1336 -  a cell affected with a starting and ending incident times
		 * 1528:668-1336;2000:1000-1500 - list of cells affected with a list of starting and ending times
		 */
		 write cells_affected;
		list<string> incidents <- cells_affected split_with ";";
		loop incident over: incidents {
			list<string> cellmom <- incident split_with ":";
			if (length(cellmom) > 1) {
				list<string> effect_se_times <- cellmom[1] split_with "-";
				effect_time <<+ int(effect_se_times[0]);
				effect_stop  <<+ int(effect_se_times[1]);
			}
			cells_aff<<+ int(cellmom[0]);
		}
		write cells_aff;
		write effect_time;
		write effect_stop;
		
		//cells_aff <- cells_affected split_with ";";
		
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
    		temp::float(read("temp")),
    		zonal_wind::float(read("zonalwind")),
    		merid_wind::float(read("meridwind")),
    		atm_press::float(read("atmpress")),
	   		co2_column::float(read("ex67")),
	   		n2_column::float(read("ex68")),
    		height::float(read("ex2")),
    		distanceFromPlanetCenter::float(read("ex1")),
			co2_atmpres_share::float(read("ex57")),
			n2_atmpres_share::float(read("ex58")),
			h20_column::float(read("ex41")),
			h20ice_column::float(read("ex43")),
			sunMarsDist::float(read("ex12")),
			solarZenithAngle::float(read("ex56")),
			o3_columnpres::float(read("ex73")),
			airHeatCap::float(read("ex8")),
			thermFluxToSpace::float(read("ex33")),
			frozenCO2::float(read("ex35")),
			emissivity::float(read("ex31"))
    		] {
    			next_step_co2_column <- co2_column;
    			int ntmp;
				int atmp; 
				albedo <- mat_albedo[id_cell-1];
				
				pCO2 <- co2_column * ga / Pa2bar; // [bar] 
				pN2 <- n2_column * ga / Pa2bar; // [bar] 
				ppH2O <- (h20_column + h20ice_column) * ga / Pa2bar; // [bar] 
				frozenCO2 <- frozenCO2 * ga / Pa2bar; /// przejście na [bar]
			

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
    }
    
    /**
     * save data every year
     */
    reflex saving_ when: cycle > 1 and (
    	cycle mod (2 * log_every_sol) = 1 // or cycle mod 24 = 1
    	) {    	
    	loop n over: cell {
    	    save [ 
    	    	n.id_cell, n.co2_column, n.nh3_column, n.ch4_column, n.cfc_column, n.ch4_column, 
    	    	n.cfc_column, n.temp, n.Ts, n.regolithCO2inv, n.energy
    	    ] to: outdir + "/year_" + int( cycle / 1336 )  + "_sol_" + (cycle / 2) mod 668 + ".csv" rewrite: false type: "csv";
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
	geometry sh <- nil;
	float temp;
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
	float o3_columnpres;
	float next_step_co2_column;
	float n2_column;
	
	float cfc_column;
	float nh3_column;
 	float ch4_column;
 	
 	float nextstep_cfc_column;
 	float nextstep_nh3_column;
 	float nextstep_ch4_column;
 	float next_step_Ts;
 	
 	float emissivity;
 	
 	float div_co2 <- 0.0;
 	float div_nh3 <- 0.0;
 	float div_ch4 <- 0.0;
 	float div_cfc <- 0.0;
 	float div_Ts <- 0.0;
 	
	float co2_atmpres_share;
	float n2_atmpres_share;
	
	float h20_column;
	float h20ice_column;
	
	float airHeatCap;
	float thermFluxToSpace;
	
	float sunMarsDist <- 1.555252;
	float distanceFromPlanetCenter;
	
	float tH2O; //  Water vapour opacity
    float frozenCO2;
    float regolithCO2inv <- regValue; // * ga / Pa2bar; // [kg*m-2] Pr (w java): [bar] = 1000 kg m-2, za: https://doi.org/10.1038/s41550-018-0529-6
    
	float S <- 1.0; //Insolation factor 
	float albedo; // albedo
	float solarZenithAngle; // solar zenith angle
	
	float pCO2 <- 0.0; 			// CO2 pressure [Pa]
	float pN2  <- 0.2e-3;		// N2 pressure [bar]
	float pCH4 <- 0.0;			// CH4 pressure [bar]
	float pNH3 <- 0.0;			// NH3 pressure [bar]
	float pCFC <- 0.0; 			// CFC pressure [bar]
	float ppH2O;				// water + ice pressure [Pa]
	
	
	float tNH3; 				// NH3 equivalent grey opacity  
	float tCFC;					// CFC equivalent grey opacity 
	float tCO2;					// CO2 equivalent grey opacity 						
	float tN2;					// N2 equivalent grey opacity 
	float tCH4; 				// CH4 equivalent grey opacity 
	
	float Ts;					// calculated mars surface temperature
	float prevTs;				// to keep the temperature from the previous step
	
	float height_diff <- 0.0;
	

	list<int> neigh <- list_with(6,0);
	list<int> azim <- list_with(6,0); 
			
	aspect base_co2{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb(255, (co2_column - 80)/2.0, 255); // border: #black;	
	}
	aspect base_Ts{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb((Ts - 125), 255, 255); // border: #black;	
	}
	aspect base_nh3{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb((nh3_column), 255, 255); // border: #black;	
	}
	aspect plain_co2{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(0, (co2_column - 80)/2.0, 0) border: #black;	
	}
	aspect plain_Ts{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb((Ts - 125), 0, 0) border: #black;	
	}
	aspect plain_nh3{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb((nh3_column), 0, 0) border: #black;	
	}
	aspect plain_energy{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(energy/3.0, 0, -energy/3.0) border: #black;	
	}
	
	
	reflex initial when: cycle = 0 {
		neighbours <- cell where (each.id_cell in neigh);
		Ts <- temp;
		prevTs <- temp;
	}
	
	reflex count_height_diff when: cycle = 0 {
		height_diff <- 0.0;
		loop n over: neighbours {
			height_diff <-  height_diff + (n.height - height);
		}
		height_diff <- height_diff * ((1-spec)/6.0 + spec/5.0);
	}
  
	
	int find_azim_hex_index( float azim_ )
	{
		if( (azim_ >= 0 and azim_ < 30) or (azim_ >= 330 and azim_ <= 360))
		{
			return 0;
		}
		loop i from: 30 to: 270 step: 60  {
			if( azim_ >= i and azim_ < i+60 )
			{
				return (i-30)/60 + 1;
			}
		}
	}	
	reflex temp_transfer when: cycle > 0 and cycle mod 2 = 1
	{
		float neigh_Ts_sum <- 0.0;
		
		loop n over: neighbours {
			neigh_Ts_sum <- neigh_Ts_sum + n.Ts * (n.height - height) * Gamma_HT;
		}
		next_step_Ts <- div_Ts + K_Ts / delta_h * (neigh_Ts_sum / (6.0 - spec) - Ts); 
	}
	
	reflex mass_transfer when: cycle > 0 and cycle mod 2 = 1 and cycle > 5 
	{
		float delta_co2_column <- 0.0;
		float neigh_co2_sum <- 0.0;
		
		loop n over: neighbours {
			neigh_co2_sum <- neigh_co2_sum + n.co2_column;

			delta_co2_column <- delta_co2_column
								- K_Co2_Temp  * 1.0/(6.0 - spec) * (n.energy - energy) / (0.5*(n.Ts + Ts) * airHeatCap) ;
			//                    	^ tu sterujemy rozrzutem temp.
		}
		co2_column <- co2_column 
						+ K_*     (((1-spec) * 4.0/6.0 + spec * 4.0/5.0) * neigh_co2_sum - 4.0 * co2_column)  
						+ delta_co2_column 
						+ div_co2;
						
		if ( co2_column < 0.0){
			co2_column <- 0.0;
		}
	}

    /**
     * w modelu zaprooponowanym przez Chrisa wszystkie fluxy są per jednostka powierzchni!
     */
	reflex _balance when: cycle > 1 and cycle mod 2 = 0 {
		prevTs <- Ts;
		prev_energy <- energy;
		
		
		float tau;
		float InsFactor <- max(0.0,sin(latitude)*sin(marsInclination*cos(sol_lon))+cos(latitude)*cos(marsInclination*cos(sol_lon)));
		float Insolation <- sunTemperature * sqrt(0.5 * sunRadius / (mat_solmars_dist[sol_lon] * AU)) * InsFactor; 
			
		float CO2MpU <- co2_column / 1.0; 				   //  the mass per unit area of the  lower atmosphere of the cell
		float RegMpU <- regolith_thickness * regolith_density; // the mass per unit area of the subsurface   	
										
		float sublimateCO2 <- 0.0; // [J]
		float resublimateCO2 <- 0.0; // [J]
		float co2_excess <- 0.0; // [kg]					
		
		div_co2 <- 0.0;

		float PsatCO2 <- 1.2264e7 * exp(-3167.8 / Ts); //[bar] za Reference this from Eq 19 in Fanale et al. (1982) Fanale, F.P., Salvail, J.R., Banerdt, W.B. and 
													// Saunders, R.S., 1982. Mars: The regolith-atmosphere-cap system and climate 
													// change. Icarus, 50(2-3), pp.381-407. 
		float PCO2 <- co2_column * ga / Pa2bar;   // [bar] aktualne ciśnienie CO2
		float co2_excess <- 0.0; // [kg]		 // nadmiar/niedobór CO2, numerycznie			
		
		//tCO2 <- 0.004 * (co2_column * ga / delta_h2)^0.4551; // Marinova et.al 2005
		tCO2 <- 0.004 * (co2_column * ga / Pa2bar)^0.4551; // Marinova et.al 2005
		
		tau <- 1.0 + tCO2;
		
		energy <- (1 - albedo) * Insolation; 
		energy <- energy + 0.75 * Insolation * emissivity * (1.0 - albedo) * tau; 
		energy <- energy - emissivity * sigma * prevTs^4;
		
		/* 
		if (cycle > 12){
			if (PCO2 > PsatCO2){ // zbyt dużo CO2 - należy go "zabrać" z atmosfery i dodać do czapy lodowej
				co2_excess <- (PCO2 - PsatCO2) / Pa2bar * ga / 100000.0;
				
				div_co2 <- -co2_excess;
				frozenCO2 <- frozenCO2 + co2_excess;
				sublimateCO2 <- co2_excess * delta_h2 * 591; // energia sublimacji [J]
				//write "freezing[" +id_cell + "]CO2 mass = " + co2_excess * delta_h2 + ", energy = " + sublimateCO2;	
			} else { // za mało CO2 - jeśli jest czapa lodowa - zabieramy trochę CO2 z czapy
				co2_excess <- (PsatCO2 - PCO2) / Pa2bar * ga / 100000.0;
				if (frozenCO2 < co2_excess){ // jeśli czapy lodowej jest zbyt mało
					co2_excess <- frozenCO2;
				} 
				div_co2 <- co2_excess;
				frozenCO2 <- frozenCO2 - co2_excess;
				resublimateCO2 <- co2_excess * delta_h2 * 591; // energia resublimacji [J]
				//write "melting[" +id_cell + "]CO2 mass = " + co2_excess * delta_h2 + ", energy = " + resublimateCO2;
			}	
		}*/
		
		
	//	energy <- energy + sublimateCO2 - resublimateCO2;
	
		if (id_cell = 333){ eneVik1 <- energy; 	}
		if (id_cell = 3090){ eneVik2 <- energy; }
		if (id_cell = 1955){ eneCell <- energy; }
		
		Ts <- energy * (delta_t / delta_h2) * 1000.0 / ( (airHeatCap * CO2MpU + regHeatCap * RegMpU))  + prevTs + next_step_Ts; 
		//             ^ tu regulujemy rozrzutem temperatury
		//Ts <- energy * 100.0 / ( (airHeatCap * CO2MpU + regHeatCap * RegMpU))  + prevTs + next_step_Ts; 
		
	}
	reflex vik1 when: id_cell = 333 {
		co2vik1 <- co2_column;
		tempVik1 <- Ts;
	}
	reflex vik2 when: id_cell = 3090 {
		co2vik2 <- co2_column;
		tempVik2 <- Ts;
	}
	reflex tcell when: id_cell = 1955 {
		co2Cell <- co2_column;
		tempCell <- Ts;
		 
	}
	/**
	 * Incident: greenhouse gas factory
	 * every sol the amount of gas increases in cell "cell_affected"
	 */
	reflex GHGfactory when: cycle > 0 and cycle mod 2 = 0 and impact_model = 1 and false {
		div_nh3 <- 0.0;
		div_cfc <- 0.0;
		div_ch4 <- 0.0;
		
		loop i from: 0 to: length(cells_aff) - 1{
			if (cycle >= effect_time[i] and cycle <= effect_stop[i] and id_cell = cells_aff[i]){ 
				div_nh3 <- nh3_const_increase / delta_h2; // add nh3_const_increase [kg] recalculated to pressure
				div_ch4 <- ch4_const_increase / delta_h2; 
				div_cfc <- cfc_const_increase / delta_h2; 				
			}
		}
	}
	reflex AsteroidImpact when: cycle > 0 and cycle mod 2 = 0 and impact_model = 2 {
		/**
		 * Ammonia might be produced on Mars biologically; 
		 * Zubrin has also proposed importing 1 - 1000 u bars of ammonia via comet impacts (Zubrin and McKay, 1997).
		 *  Here, we estimate the ammonia opacity as (Kuhn et al., 1979)
		 * 1000 [u bar] = 1000 * 10e-6 [bar] = 10e-3 [bar] = 10e-3 * 10e5 [Pa] = 10e2 [Pa]
		 */
		div_nh3 <- 0.0;
		div_cfc <- 0.0;
		div_ch4 <- 0.0;
		
		loop i from: 0 to: length(cells_aff) - 1{
			if (cycle = effect_time[i]  and id_cell = cells_aff[i]){ 
			div_nh3 <- nh3_abrupt_increase / delta_h2; // add nh3_const_increase [kg] recalculated to pressure
			div_ch4 <- ch4_abrupt_increase / delta_h2;
			div_cfc <- cfc_abrupt_increase / delta_h2;
			}
		}
	}
	
	
	
}

experiment main_experiment until: (cycle <= 100) 
{
	
	parameter "Start sollon" var: sollon_number min: 0 max: 45;
	parameter "Log every sol" var: log_every_sol;
	parameter "Regolith limit" var: regLimit min: 0 max: 1; 
	parameter "Regolith value [kgm-2] for each hex" var: regValue; 
	
	parameter "Incident model no" var: impact_model min: 0 max: 3;
	parameter "Aspect number: 1[co2], 2[temp]" var: aspect_mode min: 1 max: 3;
	
	parameter "Cells affected by incident, and timing of incident" var: cells_affected;
	
	parameter "NH3 constant increase (important only when incident model = 1) [kg]" var: nh3_const_increase;
	parameter "CH4 constant increase (important only when incident model = 1) [kg]" var: ch4_const_increase;
	parameter "CFC constant increase (important only when incident model = 1) [kg]" var: cfc_const_increase;
	
	parameter "NH3 abrupt increase (important only when incident model = 2) [kg]" var: nh3_abrupt_increase;
	parameter "CH4 abrupt increase (important only when incident model = 2) [kg]" var: ch4_abrupt_increase;
	parameter "CFC abrupt increase (important only when incident model = 2) [kg]" var: cfc_abrupt_increase;
	parameter "Eddy diffusion parameter for CO2 [m2s-1] " var: eddyDiffParam;
	parameter "Eddy diffusion parameter for NH3 [m2s-1] " var: eddyDiffParam_nh3;
	parameter "Eddy diffusion parameter for CFC [m2s-1] " var: eddyDiffParam_cfc;
	parameter "Eddy diffusion parameter for CH4 [m2s-1] " var: eddyDiffParam_ch4;
	parameter "Eddy diffusion calculated 0 [no] 1 [yes]" var: eddy_calculated min:0 max: 1;
	
	parameter "Temperature increase from warming (important only when incident model = 3 [K]" var: Tincrease;
	
	parameter "Temperature-energy coeff" var: Ts_ene_coeff; 
	parameter "Height-temp coeff" var: Gamma_HT; 
	parameter "Eddy Temperature coeff" var: eddyDiffParam_Ts;
	parameter "Eddy co2 coeff" var: eddy_co2_coeff;
	
	parameter "Folder for result" var: outdir;	
	
	output {
		/*
		display mars_co2 type: opengl camera_interaction: true ambient_light: 100 
		background: #black orthographic_projection: true rotate: 15.0 
		{
			species cell aspect: base_co2;
		}
		display mars_Ts type: opengl camera_interaction: true ambient_light: 100 
		background: #black orthographic_projection: true rotate: 15.0 
		{
			species cell aspect: base_Ts;
		}
		
		display mars_nh3 type: opengl camera_interaction: true ambient_light: 100 
		background: #black orthographic_projection: true rotate: 15.0 
		{
			species cell aspect: base_nh3;
		}
		*/
		display mars_plain_co2 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_co2;
		}
		/*
		display mars_plain_nh3 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_nh3;
		}*/
		display mars_plain_Ts type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_Ts;
		}	
		display mars_plain_ene type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_energy;
		}		
		 
		display chart_co2_vik refresh_every: 1
		{
			chart "Viking CO2, T" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "Viking 1 CO2" value: co2vik1;
				data "Viking 2 CO2" value: co2vik2;
				data "Viking 1 T" value: tempVik1;
				data "Viking 2 T" value: tempVik2;
				
				data "Selected cell T" value: tempCell;
				data "Selected cell co2" value: co2Cell;				
			}
		}
		display chart_ene_vik refresh_every: 1
		{
			chart "Viking Energy" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "Viking 1 Energy" value: eneVik1;
				data "Viking 2 Energy" value: eneVik2;
				data "Selected cell energy" value: eneCell;
				
			}
		}
			
		monitor "Mean MARS greenhouse temperature"                   name: mean_greenhouse_temp value: cell mean_of(each.Ts);
		monitor "Varianve of MARS greenhouse temperature"            name: var_greenhouse_temp value: cell variance_of(each.Ts);
		monitor "Max MARS greenhouse temperature"            		 name: max_greenhouse_temp value: cell max_of(each.Ts);
		monitor "Min MARS greenhouse temperature"            		 name: min_greenhouse_temp value: cell min_of(each.Ts);
				
		monitor "Sum of NH3" name: sumNH3 value: cell sum_of(each.nh3_column);
		monitor "Sum of CH4" name: sumCH4 value: cell sum_of(each.ch4_column);
		monitor "Sum of CFC" name: sumCFC value: cell sum_of(each.cfc_column);
		monitor "Sum of CO2" name: sumCO2 value: cell sum_of(each.co2_column);
		monitor "Number of hexes biologically habitable (T > -25 C)" name: biol_habitable_hex value: cell count(each.Ts > 248.15);
		monitor "Number of hexes with unfreezed water (T > 0 C)"     name: biosphere_hex      value: cell count(each.Ts > 273.15);

		monitor "Viking 1 CO2 pressure [kg m-2]"     name: co2v1 value: co2vik1;
		monitor "Viking 2 CO2 pressure [kg m-2]"     name: co2v2 value: co2vik2;
		monitor "Viking 1 T [K]"     name: Tv1 value: tempVik1;
		monitor "Viking 2 T [K]"     name: Tv2 value: tempVik2;
		monitor "Viking 1 Energy [J]"     name: Ev1 value: eneVik1;
		monitor "Viking 2 Energy [J]"     name: Ev2 value: eneVik2;
		
	}
}