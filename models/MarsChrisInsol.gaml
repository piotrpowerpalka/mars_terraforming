/**
* Name: MARSCAHregress
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
* 
* eddy_co2_coeff - rozrzut CO2
*/


model MARSChrisInsol

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
    float Ts_ene_coeff <- 1e5;
    float Gamma_HT <- 0.005; // Kelvin / m - potential temperature parameter 
	float ecc <- 0.09341233; // eccentricity 
	
    
    float prevTotalNH3 <- 0.0;
    float co2vik1 <- 0.0;
    float co2vik2 <- 0.0;
    float tempVik1 <- 0.0;
    float tempVik2 <- 0.0;
    float eneVik1 <- 0.0;
    float eneVik2 <- 0.0;    
    
    float co2northPole <- 0.0;
    float co2southPole <- 0.0;
    float tempNorthPole <- 0.0;
    float tempSouthPole <- 0.0;
    float co2northPoleTmp <- 0.0;
    float tempNorthPoleTmp <- 0.0;
    float co2southPoleTmp <- 0.0;
    float tempSouthPoleTmp <- 0.0;
    
    
    float co2EquCell_9 <- 0.0;   
    float co2EquCell_10 <- 0.0;
    float co2EquCell_11 <- 0.0;
    float co2EquCell_12 <- 0.0;
    float co2EquCell_637 <- 0.0;
    float co2EquCell_1065 <- 0.0;
    float co2EquCell_2442 <- 0.0;
    float co2EquCell_2455 <- 0.0;
    float co2EquCell_3459 <- 0.0;
    
    float tempEquCell_9 <- 0.0;   
    float tempEquCell_10 <- 0.0;
    float tempEquCell_11 <- 0.0;
    float tempEquCell_12 <- 0.0;
    float tempEquCell_637 <- 0.0;
    float tempEquCell_1065 <- 0.0;
    float tempEquCell_2442 <- 0.0;
    float tempEquCell_2455 <- 0.0;
    float tempEquCell_3459 <- 0.0;
    
    float tempCell0 <- 0.0;
    float tempCell15 <- 0.0;
    float tempCell30 <- 0.0;
    float tempCell45 <- 0.0;
    float tempCell60 <- 0.0;
    
    float insCell0 <- 0.0;
    float insCell15 <- 0.0;
    float insCell30 <- 0.0;
    float insCell45 <- 0.0;
    float insCell60 <- 0.0;
    
    
    float insolNorthPole <- 0.0;
    float insolSouthPole <- 0.0;
    float insolVik1 <- 0.0;
    float insolVik2 <- 0.0;
    
    float eneNorthPole <- 0.0;
    float eneSouthPole <- 0.0;
    float frozNorthPole <- 0.0;
    float frozSouthPole <- 0.0;
    
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
	
	float eddyDiffParam <- 0.05; // Eddy diffusion parameter [m2s-1] (approx.)
	 
		// dla < 0.5 praktycznie nie ma przepływu CO2 między hexami (dla ene_coeff = 1)
		
	float eddyDiffParam_cfc <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_ch4 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_nh3 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_Ts <- 0.05; // Eddy diffusion parameter [m2s-1] (approx.)  - temperature
	    // im mniej, tym bardziej temperatura rośnie, im wiecej, spada
	float eddy_co2_coeff   <- 5000;
	
	
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
	float K_Co2_Temp <- eddy_co2_coeff * delta_t / delta_h2;
	
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
    
    float regolith_thickness <- 5; // [m] grubość warstwy regolitu, która nagrzewa się od słonca (0.001m = 1mm)
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
    		height::float(read("ex2")),
			sunMarsDist::float(read("ex12")),
			airHeatCap::float(read("ex8")),
			frozenCO2::float(read("ex35")),
			emissivity::float(read("ex31"))
    		] {
    			next_step_co2_column <- co2_column;
    			int ntmp;
				int atmp; 
				albedo <- mat_albedo[id_cell-1];
				
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
	float insolat <- 0.0;
	
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
	aspect plain_insolation{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb((insolat), insolat, insolat) border: #black;	
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
	
	/**
	 * argument sin w GAMA jest podawany w stopniach
	 */
	float Insol(float sollon, float lat) {
		float S0 <- 590;		// stała słoneczna
		float e <- 0.09341233;		// mimośród (spłaszczenie orbity)
		float nachylenieOsi <- 24.936;
		float theta <- (sollon - 248) mod 360;
		float fi <- lat;
		if (fi = 90.0)  { fi <- 89.9; } 
		if (fi = -90.0) { fi <- -89.9; }
		
		float sin1 <- sin(fi);
		float kat <- nachylenieOsi * sin(sollon);
		float sin2 <- sin(nachylenieOsi) * sin(sollon);
		float cos1 <- cos(fi);
		float cos2 <- cos(kat);
		
		float zachSlonca <- -1.0 * tan(fi) * tan(kat);
		float liczbaGodzSlonecznych <- 0.0;
		
		if (zachSlonca >= 1) {
			liczbaGodzSlonecznych <- 0.0;
			zachSlonca <- 0.0;	
		}		
		else if (zachSlonca <= -1) {
			liczbaGodzSlonecznych <- 24.0;
			zachSlonca <- 180.0;
		}
		else {
			zachSlonca <- acos(zachSlonca);
			liczbaGodzSlonecznych <- zachSlonca * 2.0/15.0;
		}
		float cos3 <- cos(zachSlonca);
		float sin3 <- sin(zachSlonca);
		
		float pNawias <- S0 * (1 + e * cos(theta))^2 / (1 - e^2)^2;
		float dNawias <- sin1 * sin2 * zachSlonca * 2 * #pi / 180 + cos1 * cos2 * sin3;
		
		if (dNawias < 0 ) {
			dNawias <- 0;
		}
		float Hobh <- (24/#pi) * pNawias * dNawias;
		
		
		return Hobh/24;
	}
	
	reflex tmp when: cycle mod 2 = 0 {
		insolat <- Insol(sol_lon, latitude);
	}

    /**
     * w modelu zaprooponowanym przez Chrisa wszystkie fluxy są per jednostka powierzchni!
     */
	reflex _balance when: cycle > 1 and cycle mod 2 = 0  and false{
							
	   	float Ins <- max(0.0, Insol(sol_lon, latitude));
			
		if (id_cell = 333){ eneVik1 <- energy;  insolVik1 <- Ins;	}
		if (id_cell = 3090){ eneVik2 <- energy; insolVik2 <- Ins; }
		if (id_cell = 58)  {insolNorthPole <- Ins;  eneNorthPole <- energy; frozNorthPole <- frozenCO2; }
		if (id_cell = 1770)  {insolSouthPole <- Ins; eneSouthPole <- energy;  frozSouthPole <- frozenCO2;}
		
		
		Ts <- ((1 - albedo) * Ins / (emissivity * sigma) )^(1/4);
		//Ts <- energy * (delta_t / delta_h2) * Ts_ene_coeff / ( (airHeatCap * CO2MpU + regHeatCap * RegMpU))  + prevTs + next_step_Ts; 
		//             ^ tu regulujemy rozrzutem temperatury
	
		if (id_cell = 637) {insCell0 <- Ins; tempCell0  <- Ts; }
		if (id_cell = 597) {insCell15 <- Ins; tempCell15  <- Ts; }
		if (id_cell = 5) {insCell30 <- Ins; tempCell30  <- Ts; }
		if (id_cell = 163) {insCell45 <- Ins; tempCell45  <- Ts; }
		if (id_cell = 127) {insCell60 <- Ins; tempCell60  <- Ts; }
		
	}
	reflex vik1 when: id_cell = 333 {
		co2vik1 <- co2_column;
		tempVik1 <- Ts;
	}
	reflex vik2 when: id_cell = 3090 {
		co2vik2 <- co2_column;
		tempVik2 <- Ts;
	}
	reflex north_pole when: cycle mod 2 = 1 and (id_cell = 58) 
	 // or id_cell = 68 or 
	//	id_cell = 2850 or id_cell = 2859 or
	//	id_cell = 2860 or id_cell = 2870 or id_cell = 2871)
	 {
		//co2northPoleTmp <- co2northPoleTmp + co2_column;
		//tempNorthPoleTmp <- tempNorthPoleTmp + Ts;
		co2northPole <- co2_column;
		tempNorthPole <- Ts;		 
	}
//	reflex north_pole_0 when: cycle mod 2 = 0 {
//		co2northPole <- co2northPoleTmp / 7.0;
//		tempNorthPole <- tempNorthPoleTmp / 7.0 ;
//		co2northPoleTmp <- 0.0;
	//	tempNorthPoleTmp <- 0.0;
//		 
//	}
	reflex south_pole when: cycle mod 2 = 1 and (id_cell = 1770)
	// or id_cell = 1771  or id_cell = 1788
	//	 or id_cell = 1789 or id_cell = 1790 or id_cell = 1844 or id_cell = 1854)
	 {
		//co2southPoleTmp <- co2southPoleTmp + co2_column;
		//tempSouthPoleTmp <- tempSouthPoleTmp + Ts;
		co2southPole <- co2_column;
		tempSouthPole <- Ts;
		 
	}
//	reflex south_pole_0 when: cycle mod 2 = 0  {
//		co2southPole <- co2southPoleTmp / 7.0;
//		tempSouthPole <- tempSouthPoleTmp / 7.0 ;
//		co2southPoleTmp <- 0.0;
//		tempSouthPoleTmp <- 0.0;
//	 
//	}
	
	reflex equator_9_cell when: id_cell = 9 {
		co2EquCell_9 <- co2_column;
		tempEquCell_9 <- Ts;
	}
	reflex equator_10_cell when: id_cell = 10 {
		co2EquCell_10 <- co2_column;
		tempEquCell_10 <- Ts;
	}
	reflex equator_11_cell when: id_cell = 11 {
		co2EquCell_11 <- co2_column;
		tempEquCell_11 <- Ts;
	}
	reflex equator_12_cell when: id_cell = 12 {
		co2EquCell_12 <- co2_column;
		tempEquCell_12 <- Ts;
	}
	reflex equator_637_cell when: id_cell = 637 {
		co2EquCell_637 <- co2_column;
		tempEquCell_637 <- Ts;
	}
	reflex equator_1065_cell when: id_cell = 1065 {
		co2EquCell_1065 <- co2_column;
		tempEquCell_1065 <- Ts;
	}
	reflex equator_2442_cell when: id_cell = 2442 {
		co2EquCell_2442 <- co2_column;
		tempEquCell_2442 <- Ts;
	}
	reflex equator_2455_cell when: id_cell = 2455 {
		co2EquCell_2455 <- co2_column;
		tempEquCell_2455 <- Ts;
	}
	reflex equator_3459_cell when: id_cell = 3459 {
		co2EquCell_3459 <- co2_column;
		tempEquCell_3459 <- Ts;
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
	/* 	display mars_plain_co2 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_co2;
		}*/
		display mars_plain_insol type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_insolation;
		}
		
		
		/*
		display mars_plain_nh3 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_nh3;
		}
		display mars_plain_Ts type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_Ts;
		}	
		display mars_plain_ene type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_energy;
		}*/		
		 
		display chart_co2_poles refresh_every: 1
		{
			chart "CO2" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "North Pole CO2" value: co2northPole;
				data "South Pole CO2" value: co2southPole;
				data "Equator 9 CO2" value: co2EquCell_9;
				data "Vik1 CO2" value: co2vik1;
				data "Vik2 CO2" value: co2vik2;
				
				/*data "Equator 10 CO2" value: co2EquCell_10;
				data "Equator 11 CO2" value: co2EquCell_11;
				data "Equator 12 CO2" value: co2EquCell_12;
				data "Equator 637 CO2" value: co2EquCell_637;
				data "Equator 1065 CO2" value: co2EquCell_1065;
				data "Equator 2442 CO2" value: co2EquCell_2442;
				data "Equator 2455 CO2" value: co2EquCell_2455;
				data "Equator 3459 CO2" value: co2EquCell_3459;
				*/
			}
		}
		display chart_insol_poles refresh_every: 1
		{
			chart "Insolation" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "North Pole Insolation" value: insolNorthPole;
				data "South Pole Insolation" value: insolSouthPole;
				data "0 Insolation" value: insCell0;
				data "15 Insolation" value: insCell15;
				data "30 Insolation" value: insCell30;
				data "45 Insolation" value: insCell45;
				data "60 Insolation" value: insCell60;
			}
		}
		display chart_ene_poles refresh_every: 1
		{
			chart "Energy" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "North Pole Energy" value: eneNorthPole;
				data "South Pole Energy" value: eneSouthPole;	
				data "Vik 1 Energy" value: eneVik1;	
				data "Vik 2 Energy" value: eneVik2;	
							
			}
		}
		display chart_frozen_poles refresh_every: 1
		{
			chart "Frozen CO2" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "North Pole frozen CO2" value: frozNorthPole;
				data "South Pole frozen CO2" value: frozSouthPole;	
							
			}
		}
		display chart_temp_poles refresh_every: 1
		{
			chart "Temp" type: series size: { 1, 1.0 } position: { 0, 0 } 
			{
				data "North Pole temp" value: tempNorthPole;
				data "South Pole temp" value: tempSouthPole;
				data "0 Temp" value: tempCell0;
				data "15 Temp" value: tempCell15;
				data "30 Temp" value: tempCell30;
				data "45 Temp" value: tempCell45;
				data "60 Temp" value: tempCell60;
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
		
		/*monitor "North Pole CO2 pressure [kg m-2]"     name: co2NP value: co2northPole;
		monitor "South Pole CO2 pressure [kg m-2]"     name: co2SP value: co2southPole;
		
		monitor "Equator_9 CO2 pressure [kg m-2]"     name: co2Equ9 value: co2EquCell_9;
		monitor "Equator_10 CO2 pressure [kg m-2]"     name: co2Equ10 value: co2EquCell_10;
		monitor "Equator_11 CO2 pressure [kg m-2]"     name: co2Equ11 value: co2EquCell_11;
		monitor "Equator_12 CO2 pressure [kg m-2]"     name: co2Equ12 value: co2EquCell_12;
		monitor "Equator_637 CO2 pressure [kg m-2]"     name: co2Equ637 value: co2EquCell_637;
		monitor "Equator_1065 CO2 pressure [kg m-2]"     name: co2Equ1065 value: co2EquCell_1065;
		monitor "Equator_2442 CO2 pressure [kg m-2]"     name: co2Equ2442 value: co2EquCell_2442;
		monitor "Equator_2455 CO2 pressure [kg m-2]"     name: co2Equ2455 value: co2EquCell_2455;
		monitor "Equator_3459 CO2 pressure [kg m-2]"     name: co2Equ3459 value: co2EquCell_3459;
	
		monitor "North Pole temperature [K]"     name: tempNP value: tempNorthPole;
		monitor "South Pole temperature [K]"     name: tempSP value: tempSouthPole;
		
		monitor "Equator_9 temperature [K]"     name: tempEqu9 value: tempEquCell_9;
		monitor "Equator_10 temperature [K]"     name: tempEqu10 value: tempEquCell_10;
		monitor "Equator_11 temperature [K]"     name: tempEqu11 value: tempEquCell_11;
		monitor "Equator_12 temperature [K]"     name: tempEqu12 value: tempEquCell_12;
		monitor "Equator_637 temperature [K]"     name: tempEqu637 value: tempEquCell_637;
		monitor "Equator_1065 temperature [K]"     name: tempEqu1065 value: tempEquCell_1065;
		monitor "Equator_2442 temperature [K]"     name: tempEqu2442 value: tempEquCell_2442;
		monitor "Equator_2455 temperature [K]"     name: tempEqu2455 value: tempEquCell_2455;
		monitor "Equator_3459 temperature [K]"     name: tempEqu3459 value: tempEquCell_3459;
		*/
	}


}