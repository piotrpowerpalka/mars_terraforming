/**
* Name: SimpleDiurnal
* Based on the internal empty template. 
* Author: piotr
* Tags: 
*/


model SimpleDiurnal

/* Insert your model definition here */
global{
	int step_sol <- 24;	// number of iterations per sol
	int sol_year <- 0; // update: (cycle/(2 * step_sol)) mod martianYear; // number of sol in martian year
    int sol_lon  <- 0; // update: int(sol_year / 668.0 * 360.0); // recalculation sol no -> solar logitude
    int year <- 0; // update: int( cycle/(2 * step_sol)) - sol_year)/martianYear;   
    int hour <- 0;
    int martianYear <- 668; // martian year [sol]
	float nachylenieOsi <- 24.936;
	float S0 <- 589;		// solar constant 
	float e <- 0.09341233;		// mimośród (spłaszczenie orbity)
	float liczbaGodzSlonecznych <- 0.0;
	float wschS;
	float zachS;
	float fi;
	
	
	float temp <- 140.0;
	float energy <- 0.0; // energia w równaniu energy balance
	float prev_energy <- 0.0;
	
	float longitude;
	float latitude;
	float height;
	float Ts;
	
    float kat;
	float zachSlonca;
	float OMEGA;
	float insol;	
	
	float dustSunGlasses <- 0.0;
	float pCO2 <- 0.0; // update: atmCO2 * exp(- height / 10000); // 0.006 * exp(- height / 10000); // Pressure in [bar]
	float PsatCO2 <- 0.0; // update: 1.2264e7 * exp(-3167.8 / Ts); //[bar] za Reference this from Eq 19 in Fanale et al. (1982) Fanale, F.P., Salvail, J.R., Banerdt, W.B. and 
	
	reflex day when: cycle mod 24 = 0 {
		fi <- latitude;
		if (fi = 90.0)  { fi <- 89.99; } 
		if (fi = -90.0) { fi <- -89.99; }
		
		kat <- nachylenieOsi * sin(sol_lon);
		zachSlonca <- -1.0 * tan(fi) * tan(kat);
		
		OMEGA <- 0.0;
		
		
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
		wschS <- 12.0 - OMEGA/15.0;  // godziny 0-24
		zachS <- 12.0 + OMEGA/15.0;  // godziny 0-24
		
	}

	reflex update {
		sol_year <- (cycle/(step_sol)) mod martianYear; // number of sol in martian year
   		sol_lon  <- int(sol_year / 668.0 * 360.0); // recalculation sol no -> solar logitude
   		hour <- (cycle) mod step_sol;  
   		year <-  int( (cycle/( step_sol) - sol_year)/martianYear );
    
	}	
	
	
	reflex insolation {
		
				
		insol <- (1.0/#pi) * S0  * ((1 + e * cos((sol_lon - 248) mod 360))^2) / ((1 - e^2)^2) 
			* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sol_lon) 
					* 15.0 * 2.0 * #pi / 360.0 + cos(fi) * cos(kat) 
					* ( sin(longitude + 15.0 * (hour + 1.0)) - sin(longitude + 15.0 * hour) )
			     ) ; 		
		
	}
}

experiment main_experiment until: (cycle > 6680) 
{
	parameter "Latitude"  var: latitude min: -90 max: 90;
	parameter "Longitude" var: longitude min: -180 max: 180;
	
	
	output {
			display chart refresh: every(24#cycles) {
			chart "Slonce" type: series background: #white style: exploded {
				data "Wschod" value: wschS color: #red;
				data "Zachod" value: zachS color: #blue;
				
				}  
		}
		display chart2 refresh: every(1#cycles) {
			chart "Insolation" type: series background: #white style: exploded {
				data "Insol" value: insol color: #red;				
				}  
		}
	}
}