/**
* Name: MARSCAHregress
* Based on the internal empty template. 
* Author: Piotr Pałka
* Tags: 
* 
* eddy_co2_coeff - rozrzut CO2
*/


model MARSSingleSol


/* Insert your model definition here */

global {
	int aspect_mode <- 1;				// tryb wyświetlania
	int sollon_number <- 0;
    int model_number <- 3;
    
    int numOfHexes <- 4002;
    int logRegress <- 0;
    float Td <- 30.0; // Temperature increment to outgas 1/e of regolith
    float Tincrease <- 0.0;
    
    int GHGregress <- 0;
    float Gamma_HT <- 0.005; // Kelvin / m - potential temperature parameter 
	float ecc <- 0.09341233; // eccentricity 
	
    float eneCell <- 0.0;
    float co2Cell <- 0.0;
    float tempCell <- 0.0;
    
	string outdir <- "../results/";
	
	float Rh <- 0.7; // Relative humidity
    float Rgas <- 8.314; // Gas constant for water
    float Lheat <- 43655.0; // Latent heat (?)
    float P0 <- 1.4E6; //Reference pressure (ext(19)?)
    
	int log_every_sol <- 1;
	
	float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
    float sol_const <- 589.0; 		//Present day martian solar constant [Wm-2]
	float ga <- 3.72076; 			// Mars gravitional acceleration [ms-2] [Nkg-1]
	float AU <- 149597870700.0;		// astronomical unit [m]
	float sunTemperature <- 5780.0; 	// sun temperature [K]
	float sunRadius <- 695700000.0;	// sun radius [m]
	float marsRadius <- 3396200.0;	// mars radius [m]
	float marsArea <- 1.448e12;		// mars area [m2]\
	
	float delta_t <- 88775.0;		// 1 martian sol [s]
	float delta_h2 <-  marsArea / numOfHexes;	// area of single hex
	float delta_h <- delta_h2 / #pi;			// radius of single hex (approx.)
	
	int martianYear <- 668; // martian year [sol]
	float Pa2bar <- 1e5;			// pascal to bar
	
	// albedo
	file csv_albedo <- csv_file("../includes/heksy_albedo.csv",";",float,true);
	matrix<float> mat_albedo <- matrix<float>(csv_albedo.contents);
	
	int sol_year <- 0 update: 591;
    int sol_lon  <- 0 update: int(sol_year / 668.0 * 360.0);
    
    
	float mean_greenhouse_temp <- 0.0;  // mean greenhouse temperature
	float var_greenhouse_temp <- 0.0;
	int biol_habitable_hex <- 0;		// number of hexes biologically habitable (T > -25 C)
	int biosphere_hex <- 0; 			// number of hexes unfreezed water (T > 0 C)
	
    float marsInclination <- 25.19; // [deg] 
    
    float S0 <- 590;		// stała słoneczna
	float e <- 0.09341233;		// mimośród (spłaszczenie orbity)
	float nachylenieOsi <- 24.936;
	
	float CO2freezingCoeff <- 0.1; // percentage part of CO2 being sublimated per sol
    float CO2defreezingCoeff <- 0.05; // percentage part of CO2 being resublimated per sol
    
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
    		height::float(read("ex2")),
			sunMarsDist::float(read("ex12"))
			//frozenCO2::float(read("ex35"))
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
    reflex saving_ when: cycle > 1 and cycle mod 2 = 1{    	
    	loop n over: cell  {
    	    //save [ 
    	    //	n.id_cell, n.co2_column, // n.nh3_column, n.ch4_column, n.cfc_column, n.ch4_column, n.cfc_column, 
    	    //	n.temp, n.Ts, n.insol, n.energy
    	    //] to: outdir + "/year_" + int( cycle / 1336 )  + "_sol_" + (cycle / 2) mod 668 + ".csv" rewrite: false type: "csv";
    	    write "" + ( (cycle -1) / 2 ) + ";" + n.id_cell+ ";" + n.latitude + ";" + n.longitude +
    	         ";" +  n.Ts + ";" + n.insol + ";" + n.energy  + ";" +  n.pCO2 + ";" +  n.PsatCO2 +
    	         // + ";" + n.div_co2 +
    	         ";" + n.frozenCO2;   
    	    	
    	    
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
	
	float atm_press;
	int spec;
	
	float co2_column;
	float next_step_co2_column;
 	float next_step_Ts;
 	
 	float emissivity;
 	
 	float div_co2 <- 0.0;
 	float div_Ts <- 0.0;
 	
	float co2_atmpres_share;
	float sunMarsDist <- 1.555252;
	float distanceFromPlanetCenter;
	float frozenCO2 <- 0.0;
    
	float S <- 1.0; //Insolation factor 
	float albedo; // albedo
	float solarZenithAngle; // solar zenith angle
	
	float pCO2 <- 0.0; 			// CO2 pressure [Pa]	
	float PsatCO2 <- 0.0;
	
	
	float tCO2;					// CO2 equivalent grey opacity 						
	
	float Ts;					// calculated mars surface temperature
	float prevTs;				// to keep the temperature from the previous step
	
	float height_diff <- 0.0;
	
	float insol;
		
	float kat;
	float zachSlonca;
	float OMEGA;
	float tau;
	
	float co2_excess <- 0.0;
	float resublimateCO2 <- 0.0;
	float sublimateCO2 <- 0.0;
	
	list<int> neigh <- list_with(6,0);
	list<int> azim <- list_with(6,0); 
			
	aspect base_co2{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb(255, (co2_column - 80)/2.0, 255); // border: #black;	
	}
	aspect base_Ts{
		draw sphere(2) at: {50*x+50,50*y+50,50*z+50} color: rgb((Ts - 125), 255, 255); // border: #black;	
	}
	aspect plain_co2{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(0, (co2_column - 80)/2.0, 0) border: #black;	
	}
	aspect plain_Ts{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(Ts, Ts, Ts) border: #black;	
	}
	aspect plain_energy{
		draw circle(1.0) at: {(longitude+180.0)/3.6, (latitude+90.0)/1.8} color: rgb(energy/3.0, 0, -energy/3.0) border: #black;	
	}
	
	
	reflex initial when: cycle = 0 {
		neighbours <- cell where (each.id_cell in neigh);
		Ts <- temp;
		prevTs <- temp;
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
     * w modelu zaprooponowanym przez Chrisa wszystkie fluxy są per jednostka powierzchni!
     */
	reflex _balance when: cycle > 1 and cycle mod 2 = 0 {
							
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
		
		insol <- (24/#pi) * S0 * ((1 + e * cos((sol_lon - 248) mod 360))^2) / ((1 - e^2)^2) 
					* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sol_lon) * OMEGA * 2 * #pi / 360.0 + cos(fi) * cos(kat) * sin(OMEGA)) / 24;
		
		//Ts <- ((1 - albedo) * insol / (emissivity * sigma) )^(1/4);
		Ts <- ((1 - 0.2) * insol / (0.9 * sigma) )^(1/4);
		if (Ts = 0) { Ts <- 140.0; } // twarde ograniczenie na temp
		
		PsatCO2 <- 1.2264e7 * exp(-3167.8 / Ts); //[bar] za Reference this from Eq 19 in Fanale et al. (1982) Fanale, F.P., Salvail, J.R., Banerdt, W.B. and 
													// Saunders, R.S., 1982. Mars: The regolith-atmosphere-cap system and climate 
													// change. Icarus, 50(2-3), pp.381-407.
													
		pCO2 <- 0.006 * exp(- height / 10000); // Pressure in [bar]
		div_co2 <- 0.0;
		
		tCO2 <- 0.004 * (pCO2 / Pa2bar * ga )^0.4551; // Marinova et.al 2005 [Pa = kg*m-2 * N*kg-1 
		tau <- 1.0 + tCO2;
		
			if (pCO2 - frozenCO2 > PsatCO2){ // zbyt dużo CO2 - należy go "zabrać" z atmosfery i dodać do czapy lodowej
										// energia hexa maleje
				co2_excess <- (pCO2- frozenCO2) * CO2freezingCoeff; // [bar] 
				
				div_co2 <- -co2_excess; // [bar]
				frozenCO2 <- frozenCO2 + co2_excess; // [bar]
				//resublimateCO2 <- co2_excess * delta_h2 * 591; // energia sublimacji [J]
					
			} else if (pCO2 - frozenCO2 < PsatCO2){ // za mało CO2 - jeśli jest czapa lodowa - zabieramy trochę CO2 z czapy
				co2_excess <- frozenCO2 * CO2defreezingCoeff; // [bar]
				div_co2 <- co2_excess;
				frozenCO2 <- frozenCO2 - co2_excess; // w [bar] !!!
				//sublimateCO2 <- co2_excess * delta_h2 * 591 ; // energia resublimacji [J]				
			}
		
	}
}

experiment main_experiment until: (cycle <= 100) 
{
	
	parameter "Start sollon" var: sollon_number min: 0 max: 45;
	parameter "Log every sol" var: log_every_sol;
	
	
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

			
		monitor "Mean MARS greenhouse temperature"                   name: mean_greenhouse_temp value: cell mean_of(each.Ts);
		monitor "Varianve of MARS greenhouse temperature"            name: var_greenhouse_temp value: cell variance_of(each.Ts);
		monitor "Max MARS greenhouse temperature"            		 name: max_greenhouse_temp value: cell max_of(each.Ts);
		monitor "Min MARS greenhouse temperature"            		 name: min_greenhouse_temp value: cell min_of(each.Ts);
				
		monitor "Sum of CO2" name: sumCO2 value: cell sum_of(each.co2_column);
		monitor "Number of hexes biologically habitable (T > -25 C)" name: biol_habitable_hex value: cell count(each.Ts > 248.15);
		monitor "Number of hexes with unfreezed water (T > 0 C)"     name: biosphere_hex      value: cell count(each.Ts > 273.15);


	}


}