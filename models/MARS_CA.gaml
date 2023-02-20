/**
* Name: MARSCA
* Based on the internal empty template. 
* Author: Piotr Palka
* Tags: 
*/


model MARSCA

/* Insert your model definition here */

global {
	
	init {
		
		float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
        int Sm <- 589; 						//Present day martian solar constant
        
        create cell from: csv_file("../includes/outpkty_z_neigh_4002.csv", ";", true) with:
    	[
    		id_cell::int(read("id")),
    		x::int(read("x")),
	    	y::int(read("y")),
    		z::int(read("z")),     	
    		longitude::float(read("longitude")),
   		 	latitude::float(read("lattitude")),
    		spec::int(read("spec")),
		    
    		neigh[0]::int(read("n1")),
    		neigh[1]::int(read("n2")),
    		neigh[2]::int(read("n3")),
    		neigh[3]::int(read("n4")),
    		neigh[4]::int(read("n5")),
    		neigh[5]::int(read("n6")),
    		
    		azim[0]::int(read("a1")),
    		azim[1]::int(read("a2")),
    		azim[2]::int(read("a3")),
    		azim[3]::int(read("a4")),
    		azim[4]::int(read("a5")),
    		azim[5]::int(read("a6")),
    		
    		temp::float(read("temperture")),
    		zonal_wind::float(read("zonal_wind")),
    		meridional_wind::float(read("meridional_wind")),
    		atm_pres::float(read("atm_ressure")),
    		
    		co2_column::float(read("ex67")),
    		n2_column::float(read("ex58")),
    		
    		a::float(read("ex33"))  // ? ex(31) + ex(32) + ex(33)? - albedo
    	];
        
	}
}
 	
species cell {
	list<cell> neighbors;
	int id_cell;
	int x;
	int y;
	int z;
	geometry sh <- nil;
	float temp;
		
	float longitude;
	float latitude;
	float zonal_wind;
	float meridional_wind;
	float atm_press;
	int spec;
	
	float co2_column;
	
	float Tb; // Effective temperature of present day Mars
	float Ts; // Initialize surface temperature
	float tH2O; //  Water vapour opacity
		
	float Rh <- 0.7; // Relative humidity
    float Rgas <- 8.314; // Gas constant for water
    float Lheat <- 43655.0; // Latent heat (?)
    float P0 <- 1.4E6; //Reference pressure (ext(19)?)
    
    float regolith; // Regolith capacity of CO2
	float Pr; // Regolith inventory of CO2
	float totCO2; // Total CO2 inventory of Mars
	float Td; // Temperature increment to outgas 1/e of regolith
	float S; //Insolation factor 
	float a; // albedo
	float pCO2; //
	float pN2;
	float pCH4;
	float pNH3;
	float pCFC; // Partial pressures
    	
	float Tp;
	float Tt; // Effective, global surface, polar and tropical temperatures
	float C; // Normalization constant for regolith calculations
	float dT; // Delta T
	
	list<int> neigh; 
	//neigh <- list_with(6,0);
	
	list<int> azim; 
	//azim <- list_with(6,0);
	
	aspect base{
		if (sh != nil) {
			draw sh color: #green ;
		} 
	}
    	
    	
	reflex move_gas when: cycle > 0 {
		float wind_direction <- atan2(zonal_wind, meridional_wind);
		
//		loop idx from: 0 to: azim.length {
//			if (azim[idx] = wind_direction) {
				
//			}
//		}
	}
    	
	reflex greenhouse when: cycle > 0 {
		//tH2O <- pH2O() ^ 0.3;
		//Ts <- ( ((1 - a) * Sm * S) / (4 * sigma) ) ^ 0.25;
		//Ts <- Ts * (1 + tCO2() + tH2O() + tCH4() + tNH3() + tCFC()) ^ 0.25;
	}
	
	reflex initial when: cycle = 0 {
		if (spec = 1) {
			remove index: 6 from: neigh;
			remove index: 6 from: azim;
		}
		
		int ntmp;
		int atmp; 
		
		loop idx from: 0 to: length(neigh) { 
			loop jdx from: 0 to: length(neigh)-1 {
				if (azim[idx] > azim[jdx]) {
					atmp  <- azim[idx];
					azim[idx] <- azim[jdx];
					azim[jdx] <- atmp;
					
					ntmp  <- neigh[idx];
					neigh[idx] <- neigh[jdx];
					neigh[jdx] <- ntmp;
				}
			}
		}
		neighbors <- cell where (each.id_cell in neigh);
		
		location <- {latitude, longitude};
	}
	
	float pH2O {
		return Rh * P0 * exp(-Lheat / (Rgas * Ts));
	}
	float tCO2 {
		return 0.9 * pTot() ^ 0.45 * pCO2 ^ 0.11;
	}
	float tCH4 {
		return 0.5 * pCH4 ^ 0.278;
	}
	float tNH3 {
    	return 9.6 * pNH3 ^ 0.32;
	}
	float tCFC {
		return 1.1 ^ pCFC / (0.015 + pCFC);
	}

	float pTot {
		return pCO2 + pN2 + pCH4 + pNH3 + (pCFC / 1E5) + pH2O();
	}
    

}
experiment main_experiment until: (cycle <= 8065)
{
	output {
		display mars type: opengl ambient_light: 100
		{
			species cell aspect: base;
		}
	}
}