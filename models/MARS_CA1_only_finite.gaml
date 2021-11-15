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
    int model_number <- 0;
	string dir_res <- "finite_full";
    int height_diff_include <- 0;
    
    float b0 -> -0.04051304926747448;
	float b1 -> 1.59819953e-02;
	float b2 -> -1.20888544e-02; 
	float b3 -> -6.26106713e-05;
	float b4 -> -7.87344857e-05; 
	float b5 -> 9.47135571e-05;
	
	float a -> 0.0071;
	float b -> -0.06;
	
	float delta_h2 <- 1.448e9 / 4;
	float delta_h <- delta_h2 / 3.14;
	float K_ <- 0.2452348; // K * delta_t / delta_h2
        
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
	   		n2_column::float(read("ex58")),
    		albedo::float(read("ex33")),
    		height::float(read("ex2"))
    		] {
    			next_step_co2_column <- co2_column;
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
    
    reflex saving_ when: cycle > 0 and cycle mod 20 = 0 {    	
    	loop n over: cell {
    	    save [ 
    	    	n.id_cell, n.co2_column, n.temp, n.Ts
    	    ] to: "../results/" + dir_res + "/sol" + sollon_number + "/save_csv_" + (cycle / 2) + ".csv" rewrite: false type: "csv";
   	 	}	
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
	
	reflex move_gas when: cycle mod 2 = 1 and model_number = 0 {
		if spec = 1
		{
			next_step_co2_column <- co2_column;
			return;
		}
		
		//if height_diff_include = 1 {
		//	do height_diff_co2_column;	
		//}
		
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
	
	reflex move_gas_regression when: cycle > 0 and cycle mod 2 = 1 and model_number = 1
	{
		float dw <- atan2( zonal_wind, merid_wind ); // * 180 / 3.14;
		float xs <- sqrt( zonal_wind^2 + merid_wind^2 );
		
		//if height_diff_include = 1 {
		//	do height_diff_co2_column;	
		//}
		
		loop i from: 0 to: length(neighbours)-1
		{
			float xa <- min( (azim[i] - dw ) mod 360, 360 - (azim[i] - dw) mod 360 );
			
			float delta_co2_column <- 1/6*(b0 + b1*xa + b2*xs + b3*xa*xa + b4*xa*xs + b5*xs*xs);
			
			next_step_co2_column <- next_step_co2_column - delta_co2_column;
			
			ask neighbours[i] {
				next_step_co2_column <- next_step_co2_column + delta_co2_column;
			}
		}
	}
    		
	reflex height_diff_co2_column when: cycle > 0 and cycle mod 2 = 1 and model_number = 2
	{
		loop n over: neighbours
		{
			float height_diff <- n.height - height;
			float delta_co2 <- 1/6 * (b + a * height_diff);
			
			next_step_co2_column <- next_step_co2_column - delta_co2;
			n.next_step_co2_column <- n.next_step_co2_column + delta_co2;
		}
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
		
		float param <- (1-spec) * 2/3 + spec * 4/5; // jezeli spec to 4/5 inaczej 2/3

		if( alfa = beta ){
			next_step_co2_column <- co2_column;
		}else{
			float alfa_dw <- zonal_wind * (sin(beta) - cos(beta)) / sin(beta - alfa);
			float beta_dw <- merid_wind * (cos(alfa) - sin(alfa)) / sin(beta - alfa);
		
			float neigh_co2_sum <- 0;
			float height_diff <- 0.0;
			
			loop n over: neighbours {
				neigh_co2_sum <- neigh_co2_sum + n.co2_column;
				height_diff <-  height_diff + (n.height - height);
			}
			
			height_diff <- height_diff * ((1-spec)/6.0 + spec/5.0);
		
			next_step_co2_column <- -(beta_dw / delta_h)*(neighbours[beta_az].co2_column - co2_column) - (alfa_dw / delta_h)*(neighbours[alfa_az].co2_column - co2_column) + K_* (param*neigh_co2_sum - 4*co2_column) + co2_column;
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
	
//	reflex greenhouse when: cycle > 0 and cycle mod 2 = 0 {
//		co2_column <- next_step_co2_column;
//		
//		C <- regolith * 0.006^(-0.275) * exp(149 / Td);
//		
//		totCO2 <- Pr + pole + pCO2;
//		
//		Tb <- ((1.0 - 0.2)*(sol_const)/(4.0*sigma)) ^ 0.25;
//		Ts <- Tb;
//		Tp <- Ts - 75;
//		
//		loop lo from:1 to: 100 {
//			tH2O <- pH2O() ^ 0.3;
//			Ts <- ( ((1.0 - albedo/100.0)*(S*sol_const)/(4.0*sigma)) ^ 0.25);
//			Ts <-  Ts * (1 + tCO2() + tH2O + tCH4() + tNH3() + tCFC()) ^ 0.25;
//			
//			Tp <- Ts - 75 / (1 + 5 * pTot());
//			pCO2 <- pressCO2(); 
//		}
//		Tt <- Ts * 1.1;
//		dT <- Ts - Tb;		
//		
//		ppH2O <- pH2O();
//	}

	reflex update when: cycle > 0 and cycle mod 2 = 0 {
		co2_column <- next_step_co2_column;
	}

	reflex greenhouse when: cycle < 0 {
		Ts <- ((1.0 - 0.2)*(sol_const)/(4.0*sigma)) ^ 0.25;
		
		tH2O <- pH2O() ^ 0.3;
		
		Ts <- ( ((1.0 - albedo/100.0)*(S*sol_const)/(4.0*sigma)) ^ 0.25);
		Ts <-  Ts * (1 + tCO2() + tH2O + tCH4() + tNH3() + tCFC()) ^ 0.25;
		
		pCO2 <- pressCO2();		
	}

	float pressCO2 {
		float co2_V <- delta_h2 * co2_column;
		float co2_number_of_moles <- ( atm_press * co2_V ) / ( (Ts+273) * 8.3145 ); // temp albo Ts
		float co2_mass <- 0.001 * 44.009 * co2_number_of_moles; // kg
		return co2_mass / delta_h2;
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
	parameter "Sollon startowy" var: sollon_number min: 0 max: 45;
	parameter "Numer modelu" var: model_number min: 3;
	parameter "Uwzględnienie różnicy wysokości" var: height_diff_include min: 0 max: 1;
	parameter "Katalog z wynikami" var: dir_res;
	
	output {
		display mars type: opengl ambient_light: 100
		{
			species cell aspect: base;
		}
	}
}
