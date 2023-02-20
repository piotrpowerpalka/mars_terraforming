/**
* Name: MARSCA1
* Based on the internal empty template. 
* Author: Piotr Palka
* Tags: 
*/


model MARSCA1

/* Insert your model definition here */

global {
	
	float sigma <- 0.00000005670374419; // Stefan-Boltzmann constant
    float sol_const <- 589.0; 						//Present day martian solar constant
        float moving_percentage_of_gas <- 0.01;
	float opening_angle <- 60.0;
	string sollon_number <- 0;
	init {
		create cell from: csv_file("../includes/out_grid_4002_0h_" + sollon_number + "_sollon.csv", ";", true) with:
    	[   id_cell::int(read("id")),x::int(read("x")),y::int(read("y")),z::int(read("z")),     	
    		longitude::float(read("longitude")),latitude::float(read("lattitude")),spec::int(read("spec")),
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
	   		n2_column::float(read("ex58")),
    		albedo::float(read("ex33"))
    		] {
    			albedo <- 20.0;
				int ntmp;
				int atmp; 

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
    
    reflex step when: cycle > 0 { 
    	ask cell parallel: true{
    		do next_step;
    	}
    }
//    
//    reflex saving when: cycle = 2 {
//		save species_of(cell) to: "save_csvfile123.csv" type: "csv" header: false attributes: ["x", "y", "co2_column"];
//
//    }
    
    reflex saving_ {    	
    	loop n over: cell {
    	    save [ n.id_cell, n.x, n.y, n.z, n.n1, n.n2, n.n3, n.n4, n.n5, n.n6, n.a1, n.a2, n.a3, n.a4, n.a5, n.a6, n.temp, n.longitude, n.latitude, n.zonal_wind, n.merid_wind, n.atm_press, n.spec, n.co2_column, n.next_step_co2_column, n.n2_column, n.Tb, n.tH2O, n.Rh, n.Rgas, n.Lheat, n.P0, n.regolith, n.pole, n.Pr, n.totCO2, n.Td, n.S, n.albedo, n.pCO2, n.pN2, n.pCH4, n.pNH3, n.pCFC, n.ppH2O, n.Tp, n.Tt, n.Ts, n.C, n.dT] to: "../results/sol" + sollon_number + "/save_csv_" + cycle + ".csv" rewrite: false type: "csv";
   	 	}	
   	 }
    
}

    
species cell {
	list<cell> neighbours;
	
	int id_cell;
	int x;
	int y;
	int z;
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
	float zonal_wind;
	float merid_wind;
	float atm_press;
	int spec;
	
	float co2_column;
	float next_step_co2_column;
	float n2_column;
	
	
	float Tb; // Effective temperature of present day Mars
	float tH2O; //  Water vapour opacity
	
	float Rh <- 0.7; // Relative humidity
    float Rgas <- 8.314; // Gas constant for water
    float Lheat <- 43655.0; // Latent heat (?)
    float P0 <- 1.4E6; //Reference pressure (ext(19)?)
    
    float regolith <- 300 / 1e3; // Regolith capacity of CO2
    float pole <- 50 / 1e3;
	float Pr <- 300 / 1e3; // Regolith inventory of CO2
	float totCO2; // Total CO2 inventory of Mars
	float Td <- 30.0; // Temperature increment to outgas 1/e of regolith
	float S <- 1.0; //Insolation factor 
	float albedo; // albedo
	
	float pCO2 <- 0.0; //
	float pN2  <- 0.2 / 1e3;
	float pCH4 <- 0.0;
	float pNH3 <- 0.0;
	float pCFC <- 0.0; // Partial pressures
	float ppH2O;
	
	float Tp;
	float Tt; // Effective, global surface, polar and tropical temperatures
	float Ts;
	float C; // Normalization constant for regolith calculations
	float dT; // Delta T
	
	list<int> neigh <- list_with(6,0);
	list<int> azim <- list_with(6,0); 
			
	aspect base{
		draw circle(1) color: rgb(temp, 0,0) border: #black;
	}
	
	action next_step {
		do move_gas;
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
	
	action move_gas {
		if spec = 1
		{
			next_step_co2_column <- co2_column;
			return;
		}
		float wind_direction <- atan2(zonal_wind, merid_wind);
		
		float x1 <- wind_direction - opening_angle/2;
		
		if( x1 < 0 ){
			x1 <- 360.0 + x1;
		}
		
		float x2 <- wind_direction + opening_angle/2;
		
		if( x2 > 360.0 ){
			x2 <- x2 - 360.0;
		}
		
		int i1 <- find_azim_hex_index( x1 );
		int i2 <- find_azim_hex_index( x2 ); 
		
		float co2_to_be_moved <- moving_percentage_of_gas * co2_column;
		
		next_step_co2_column <- co2_column - co2_to_be_moved;
		
		if( i1 = i2 )
		{
			// cały gaz idzie do jednego hexa
			ask neighbours[i1] {
				next_step_co2_column <- next_step_co2_column + co2_to_be_moved;
			}
		}
		else 
		{
			ask neighbours[i1] {
				float proportion <- 0.0;
				
				if( x1 > 330 )
				{
					float proportion <- ( 360 - x1 + 30 ) / opening_angle;
				}
				else
				{
					float proportion <- abs(( 30 + i1 * 60 )-x1) / opening_angle;	
				}
				next_step_co2_column <- next_step_co2_column + proportion*co2_to_be_moved;
			}
			
			ask neighbours[i2] {
				float proportion <- 0.0;
				
				if( x2 < 30 )
				{
					float proportion <- ( x2+30 ) / opening_angle;
				}
				else
				{
					float proportion <- ( x2-30 mod 60 ) / opening_angle;	
				}
				next_step_co2_column <- next_step_co2_column + proportion*co2_to_be_moved;
			}
			
			if ( i1 < i2 )
			{
				// kierunek północny nie znajduje się między x1 i x2
				loop i from: i1+1 to: i2-1 {
					// hexy pomiędzy i1 i i2 otrzymują część gazu proporcjonalną do (60 / opening_angle) 
					ask neighbours[i] {
						 next_step_co2_column <- next_step_co2_column + (60/opening_angle)*co2_to_be_moved;
					}
				}	
			}else{
				// kierunek północny znajduje się miedzy x1 i x2
				loop i from: 0 to: length(neighbours)-1 {
					// wszystkie hexy poza tymi, znajdującymi sie miedzy i1 i i2 
					// otrzymują część gazu proporcjonalną do 60/opening_angle
					if( i < i2 or i > i1 ) {
						ask neighbours[i] {
							next_step_co2_column <- next_step_co2_column + (60/opening_angle)*co2_to_be_moved;
						}
					}
				}
			}
			
		}
	}
    	
	reflex greenhouse when: cycle > 0 {
		co2_column <- next_step_co2_column;
		
		C <- regolith * 0.006^(-0.275) * exp(149 / Td);
		
		//
		totCO2 <- Pr + pole + pCO2;
		
		Tb <- ((1.0 - 0.2)*(sol_const)/(4.0*sigma)) ^ 0.25;
		Ts <- Tb;
		Tp <- Ts - 75;
		
		loop lo from:1 to: 100 {
			tH2O <- pH2O() ^ 0.3;
			Ts <- ( ((1.0 - albedo/100.0)*(S*sol_const)/(4.0*sigma)) ^ 0.25);
			Ts <-  Ts * (1 + tCO2() + tH2O + tCH4() + tNH3() + tCFC()) ^ 0.25;
			
			Tp <- Ts - 75 / (1 + 5 * pTot());
			pCO2 <- pressCO2(); 
		}
		Tt <- Ts * 1.1;
		dT <- Ts - Tb;		
		
		ppH2O <- pH2O();
		
	}
	
	float pressCO2 {
        // Calculates pressure of CO2 after partitioning between cap
        // atmosphere and regolith
        float Pv; // CO2 polar vapour pressure
        float Pa; // Temporary variable for CO2 pressure estimate
        float X;
        float Y; // Working variables
        float top;
        float bottom; // Bisection method variables
        
        Pa <- pCO2;
        Pv <- 1.23E7 * exp(-3168 / Tp);
        
        if (Pv > Pa and pole > 0 and Pv < Pa + pole) {
        	pole <- pole - (Pv - Pa);
            Pa <- Pv;
        }
        
        if (Pv > Pa + pole and pole > 0) {
            Pa <- Pa + pole;
            pole <- 0.0;
        }
        
        if (Pv < Pa)
        {
            pole <- pole - (Pv - Pa);
            Pa <- Pv;
        }
  
        X <- Pa + Pr;
        if (X > totCO2){
        	X <- totCO2;
        }
                
        Y <- C * exp(-Tp / Td);
        top <- totCO2;
        bottom <- 0.0;
        Pa <- 0.5 * regolith;
        
        // Calculation of Pr by bisection method
        loop lop from: 1 to: 50 {
           if (Y * (Pa ^ 0.275) + Pa < X){
           		bottom <- Pa;
           } else {
           		top <- Pa;
           }
           Pa <- bottom + (top - bottom) / 2;
        }
        
         Pr <- Y * (Pa ^ 0.275);
         if (Pr > regolith) {
           		Pa <- Pa + Pr - regolith;
                Pr <- regolith;
         }  
        return Pa;
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
	
	/* Ok */
	float pH2O {
		return Rh * P0 * exp(-Lheat / (Rgas * Ts));
	}
	/* Ok */
	float tCO2 {
		return 0.9 * pTot() ^ 0.45 * pCO2 ^ 0.11;
	}
	/* Ok */
	float tCH4 {
		return 0.5 * pCH4 ^ 0.278;
	}
	/* Ok */
	float tNH3 {
    	return 9.6 * pNH3 ^ 0.32;
	}
	/* Ok */
	float tCFC {
		return 1.1 * pCFC / (0.015 + pCFC);
	}

	/* Ok */
	float pTot {
		return pCO2 + pN2 + pCH4 + pNH3 + (pCFC / 1E5) + pH2O();
	}
    	
}

experiment main_experiment until: (cycle <= 100)
{
	parameter "Procent gazu wywiewanego w pojedynczym cylku z hexa" var: moving_percentage_of_gas  min: 0.0 max: 1.0;
	parameter "Rozwarcie stożka wiatru" var: opening_angle min: 0.0 max: 360.0;
	parameter "Sollon startowy" var: sollon_number min: 0 max: 100;
	
	output {
		display mars type: opengl ambient_light: 100
		{
			species cell aspect: base;
		}
	}
}
//experiment batch_experiment type: batch repeat: 2 keep_seed: true until: ( cycle > 1000 ) {
//	parameter "Procent gazu wywiewanego w pojedynczym cylku z hexa" var: moving_percentage_of_gas  min: 0.0 max: 1.0 step: 0.25;
//	
//}
