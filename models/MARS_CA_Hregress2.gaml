/**
* Name: MARSCAHregress
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
*/


model MARSCAHregress

/* Insert your model definition here */

global {
	int aspect_mode <- 1;				// tryb wyświetlania
	int green_mode <- 5;
	int impact_model <- 1;
    float moving_percentage_of_gas <- 0.01;
    float opening_angle <- 60.0;
	int sollon_number <- 0;
    int model_number <- 0;
    int height_diff_include <- 1;
    int numOfHexes <- 4002;
    
    int cell_affected <- 1000; 				// cell affected by the effect
    int effect_time <- 30;					// cycle of the efect happened
    
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
	
	
	float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
    float sol_const <- 589.0; 		//Present day martian solar constant [Wm-2]
	float ga <- 3.72076; 			// Mars gravitional acceleration [ms-2] [Nkg-1]
	float AU <- 149597870700.0;		// astronomical unit [m]
	float sunTemperature <- 5780.0; 	// sun temperature [K]
	float sunRadius <- 695700000.0;	// sun radius [m]
	float marsRadius <- 3396200.0;	// mars radius [m]
	float marsArea <- 1.448e12;		// mars area [m2]\
	
	float eddyDiffParam <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_cfc <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_ch4 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	float eddyDiffParam_nh3 <- 1000.0;	// Eddy diffusion parameter [m2s-1] (approx.)
	
	float delta_t <- 88775.0;		// 1 martian sol [s]
	
	float delta_h2 <-  marsArea / numOfHexes;	// area of single hex
	float delta_h <- delta_h2 / #pi;			// radius of single hex (approx.)
	
	float K_ <- eddyDiffParam * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_CFC_ <- eddyDiffParam_cfc * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_CH4_ <- eddyDiffParam_ch4 * delta_t / delta_h2; // K * delta_t / delta_h2
	float K_NH3_ <- eddyDiffParam_nh3 * delta_t / delta_h2; // K * delta_t / delta_h2
	
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
			thermFluxToSpace::float(read("ex33"))
    		] {
    			next_step_co2_column <- co2_column;
    			int ntmp;
				int atmp; 
				albedo <- mat_albedo[id_cell-1];
				
				pCO2 <- co2_column * ga / Pa2bar; // [bar] 
				pN2 <- n2_column * ga / Pa2bar; // [bar] 
				ppH2O <- (h20_column + h20ice_column) * ga / Pa2bar; // [bar] 
			

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
    	cycle mod 1336 = 0 // or (cycle / 2) mod 167 = 0
    	) {    	
    	loop n over: cell {
    	    save [ 
    	    	n.id_cell, n.co2_column, n.nh3_column, n.ch4_column, n.cfc_column, n.temp, n.Ts
    	    ] to: outdir + "/year_" + int( cycle / 1336 )  + "_sol_" + (cycle / 2) mod 668 + ".csv" rewrite: false type: "csv";
   		}
    }
    
    reflex when: cycle > 1 and cycle mod 2 = 1 and cycle <= 1336 {
    	int sollon <- int((cycle-1)/1336.0 * 360.0) mod 360;
    	int solno <- int((cycle-1)/2) mod 668;
    	
    	//write "sollon = " + sollon + ", sol = " + solno;
    	loop n over: cell {
			instances[0, n.id_cell - 1] <- n.next_step_co2_column - mat_co2[n.id_cell - 1, sollon]; // TODO: dopasować do 1336
			instances[1, n.id_cell - 1] <- n.height_diff;
		}    	
		co2_fct_mat[solno] <- build(instances);
		
		//write "learnt function: " + co2_regress_fct;
		
    }
}

species cell {
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
 	
	
	float co2_atmpres_share;
	float n2_atmpres_share;
	
	float h20_column;
	float h20ice_column;
	
	float airHeatCap;
	float thermFluxToSpace;
	
	float sunMarsDist <- 1.555252;
	float distanceFromPlanetCenter;
	
	float tH2O; //  Water vapour opacity
    
    float Td <- 30.0; // Temperature increment to outgas 1/e of regolith
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
	
	
	reflex initial when: cycle = 0 {
		neighbours <- cell where (each.id_cell in neigh);
		if (spec != 1) {
    		sh <- polygon([{(neighbours[0].longitude + longitude)/2, (neighbours[0].latitude + latitude)/2}, 
				{(neighbours[1].longitude + longitude)/2, (neighbours[1].latitude + latitude)/2},
				{(neighbours[2].longitude + longitude)/2, (neighbours[2].latitude + latitude)/2},
				{(neighbours[3].longitude + longitude)/2, (neighbours[3].latitude + latitude)/2},
				{(neighbours[4].longitude + longitude)/2, (neighbours[4].latitude + latitude)/2},
				{(neighbours[5].longitude + longitude)/2, (neighbours[5].latitude + latitude)/2}
			]);
    		
    	} else {
    		sh <- polygon([{(neighbours[0].longitude + longitude)/2, (neighbours[0].latitude + latitude)/2}, 
				{(neighbours[1].longitude + longitude)/2, (neighbours[1].latitude + latitude)/2},
				{(neighbours[2].longitude + longitude)/2, (neighbours[2].latitude + latitude)/2},
				{(neighbours[3].longitude + longitude)/2, (neighbours[3].latitude + latitude)/2},
				{(neighbours[4].longitude + longitude)/2, (neighbours[4].latitude + latitude)/2}
			]);
    	}
	}
	
	reflex count_height_diff when: cycle = 0 {
		height_diff <- 0.0;
		loop n over: neighbours {
			height_diff <-  height_diff + (n.height - height);
		}
		height_diff <- height_diff * ((1-spec)/6.0 + spec/5.0);
	}
  
	reflex finite_element when: cycle > 0 and cycle mod 2 = 1 and model_number = 3
	{
		float dw <- atan2( zonal_wind, merid_wind ); // * 180 / 3.14;
		int beta_az <- length(azim)-1;
		loop az_index from: 0 to: length(azim) - 2 
		{
			if( azim[az_index] < dw and azim[az_index + 1] >  dw )
			{
				beta_az <- az_index;
			}
		}
		
		int alfa_az <- (beta_az+1) mod (length(azim)-1);
		
		float alfa <- azim[alfa_az];
		float beta  <- azim[beta_az];
		
		float param <- (1-spec) * 4.0/6.0 + spec * 4.0/5.0; // jezeli spec to 4/5 inaczej 2/3

		if( alfa = beta ){
			next_step_co2_column <- co2_column;
			nextstep_cfc_column <- cfc_column;
			nextstep_ch4_column <- ch4_column;
			nextstep_nh3_column <- nh3_column;
		} else {
			float alfa_dw <- zonal_wind * (sin(beta) - cos(beta)) / sin(beta - alfa);
			float beta_dw <- merid_wind * (cos(alfa) - sin(alfa)) / sin(beta - alfa);
		
			float neigh_co2_sum <- 0;
			float neigh_cfc_sum <- 0;
			float neigh_nh3_sum <- 0;
			float neigh_ch4_sum <- 0; 
			
			loop n over: neighbours {
				neigh_co2_sum <- neigh_co2_sum + n.co2_column;
				neigh_cfc_sum <- neigh_cfc_sum + n.cfc_column;
				neigh_ch4_sum <- neigh_ch4_sum + n.ch4_column;
				neigh_nh3_sum <- neigh_nh3_sum + n.nh3_column;
			}
			next_step_co2_column <- -(beta_dw / delta_h)*(neighbours[beta_az].co2_column - co2_column) - (alfa_dw / delta_h)*(neighbours[alfa_az].co2_column - co2_column) + K_* (param*neigh_co2_sum - 4*co2_column) + co2_column;
			nextstep_cfc_column <- -(beta_dw / delta_h)*(neighbours[beta_az].cfc_column - cfc_column) - (alfa_dw / delta_h)*(neighbours[alfa_az].cfc_column - cfc_column) + K_CFC_* (param*neigh_cfc_sum - 4*cfc_column) + cfc_column;
			nextstep_ch4_column <- -(beta_dw / delta_h)*(neighbours[beta_az].ch4_column - ch4_column) - (alfa_dw / delta_h)*(neighbours[alfa_az].ch4_column - ch4_column) + K_CH4_* (param*neigh_ch4_sum - 4*ch4_column) + ch4_column;
			nextstep_nh3_column <- -(beta_dw / delta_h)*(neighbours[beta_az].nh3_column - nh3_column) - (alfa_dw / delta_h)*(neighbours[alfa_az].nh3_column - nh3_column) + K_NH3_* (param*neigh_nh3_sum - 4*nh3_column) + nh3_column;
			
			if (height_diff_include = 1 and co2_fct_mat[(cycle/2) mod 668] != nil) {
				next_step_co2_column <- next_step_co2_column -  predict(co2_fct_mat[(cycle/2) mod 668], [height_diff]);	
			}
			nextstep_nh3_column  <- nextstep_nh3_column  -  predict(co2_fct_mat[(cycle/2) mod 668], [height_diff]);
			if (nextstep_nh3_column < 0) { nextstep_nh3_column <- 0.0; }	

			nextstep_ch4_column  <- nextstep_ch4_column  - predict(co2_fct_mat[(cycle/2) mod 668], [height_diff]);
			if (nextstep_ch4_column < 0) { nextstep_ch4_column <- 0.0; }	

			nextstep_cfc_column  <- nextstep_cfc_column  - predict(co2_fct_mat[(cycle/2) mod 668], [height_diff]);	
			if (nextstep_cfc_column < 0) { nextstep_cfc_column <- 0.0; }	
			
		}
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
	
	reflex update when: cycle > 0 and cycle mod 2 = 0 {
		co2_column <- next_step_co2_column;
		cfc_column <- nextstep_cfc_column;
		ch4_column <- nextstep_ch4_column;
		nh3_column <- nextstep_nh3_column;	
	}

	
/**
	 * model 2021-09-25 v5
	 */
	reflex greenhouse_5 when: cycle > 1 and cycle mod 2 = 0 and green_mode = 5 {
		prevTs <- Ts;
		Ts <- sunTemperature * sqrt(0.5 * sunRadius / (mat_solmars_dist[sol_lon] * AU)) * (1.0 - albedo);
		//Ts <- sunTemperature * sqrt(0.5 * sunRadius / (sunMarsDist * AU)) * (1.0 - albedo);
		
		pCO2  <- co2_column * ga / Pa2bar; // [bar] 	
		ppH2O <- (h20_column + h20ice_column) * ga / Pa2bar;
		pCH4  <- ch4_column * ga / Pa2bar; // [bar]
		pCFC  <- cfc_column * ga / Pa2bar; // [bar]
		pNH3  <- nh3_column * ga / Pa2bar; // [bar]
		
		tH2O <- ppH2O ^ 0.3;
		tCO2 <- 0.9 * (atm_press / Pa2bar) ^ 0.45 * pCO2 ^ 0.11;
		tCH4 <- 0.5 * (ch4_column / Pa2bar) ^ 0.278;
		tNH3 <- 9.6 * (nh3_column / Pa2bar) ^ 0.32;
		tCFC <- 1.1 * (cfc_column / Pa2bar) / (0.015 + (pCFC / Pa2bar) );
		
		Ts <-  Ts * (1 + tCO2+ tH2O + tCH4 + tNH3 + tCFC) ^ 0.25;
		Ts <- Ts + (prevTs * airHeatCap) * 0.000742991 - 96.96617837;
		Ts <- Ts + albedo * 100.0 * 2.062852087 - 41.02245772;
		Ts <- Ts + thermFluxToSpace * 0.209982289 - 13.16459818;
	}
	
	/**
	 * Incident: greenhouse gas factory
	 * every sol the amount of gas increases in cell "cell_affected"
	 */
	reflex GHGfactory when: cycle > 1 and cycle mod 2 = 0 and impact_model = 1 {
		if (id_cell = cell_affected){
			nh3_column <- nh3_column + nh3_const_increase / delta_h2; // add nh3_const_increase [kg] recalculated to pressure
			ch4_column <- ch4_column + ch4_const_increase / delta_h2; 
			cfc_column <- cfc_column + cfc_const_increase / delta_h2; 
		}
	}
	reflex AsteroidImpact when: cycle =  effect_time and impact_model = 2 {
		/**
		 * Ammonia might be produced on Mars biologically; 
		 * Zubrin has also proposed importing 1 - 1000 u bars of ammonia via comet impacts (Zubrin and McKay, 1997).
		 *  Here, we estimate the ammonia opacity as (Kuhn et al., 1979)
		 * 1000 [u bar] = 1000 * 10e-6 [bar] = 10e-3 [bar] = 10e-3 * 10e5 [Pa] = 10e2 [Pa]
		 */
		 if (id_cell = cell_affected){
			nh3_column <- nh3_column + nh3_abrupt_increase / delta_h2; // add nh3_const_increase [kg] recalculated to pressure
			ch4_column <- ch4_column + ch4_abrupt_increase / delta_h2;
			cfc_column <- cfc_column + cfc_abrupt_increase / delta_h2;
		}
	}
	
	reflex PoleMelting when: cycle > 1 and cycle mod 2 = 0 and impact_model = 3 {
		
	}
}

experiment main_experiment until: (cycle <= 100)
{
	parameter "Start sollon" var: sollon_number min: 0 max: 45;
	parameter "Gas transport model" var: model_number min: 3;
	parameter "Height regression include" var: height_diff_include min: 0 max: 1;
	parameter "Greenhouse model" var: green_mode min: 0 max: 5;
	parameter "Incident model no" var: impact_model min: 0 max: 3;
	parameter "Aspect number: 1[co2], 2[temp]" var: aspect_mode min: 1 max: 3;
	
	parameter "Cell affected by incident (important only when incident model = 2)" var: cell_affected;
	parameter "Cycle no of the incident (important only when incident model = 2)" var: effect_time;
	
	
	parameter "NH3 constant increase (important only when incident model = 1)" var: nh3_const_increase;
	parameter "CH4 constant increase (important only when incident model = 1)" var: ch4_const_increase;
	parameter "CFC constant increase (important only when incident model = 1)" var: cfc_const_increase;
	
	parameter "NH3 abrupt increase (important only when incident model = 2)" var: nh3_abrupt_increase;
	parameter "CH4 abrupt increase (important only when incident model = 2" var: ch4_abrupt_increase;
	parameter "CFC abrupt increase (important only when incident model = 2" var: cfc_abrupt_increase;
	
	parameter "Folder for result" var: outdir;
	
	
	output {
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
		display mars_plain_co2 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_co2;
		}
		display mars_plain_nh3 type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_nh3;
		}
		display mars_plain_Ts type: opengl ambient_light: 100 background: #white  
		{
			species cell aspect: plain_Ts;
		}		
	}
}