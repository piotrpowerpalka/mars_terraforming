/**
* Name: balanceModel
* Based on the internal empty template. 
* Author: piotr
* Tags: 
*/


model balanceModel

/* Insert your model definition here */

species cell parallel: true {
	float energy <- 0.0; // energia w równaniu energy balance
	float heat_flux;
	float Ts;
	float prevTs;
	float latitude;
	float kat;
	float nachylenieOsi;
	float sol_lon;
	float zachSlonca;
	float OMEGA;
	float S0 <- 589;
	float e;
	float pCO2;
	float pSatCO2;
	
	
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
		
		insol <- (1.0/#pi) * S0 * ((1 + e * cos((sol_lon - 248) mod 360))^2) / ((1 - e^2)^2) 
					* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sol_lon) 
					* OMEGA * 2 * #pi / 360.0 + cos(fi) * cos(kat) * sin(OMEGA)
					) ; 		
		
		prevTs <- Ts;
		bool pCO2greater <- false;
		bool check <- true;
		float co2_excess_factor <- 1.0;
		
		loop times: 10 {
			
			pCO2 <- atmCO2 * exp(- height / 10000); // 0.006 * exp(- height / 10000); // Pressure in [bar];
			if (pCO2 < 1e-8) {pCO2 <- 1e-8; }
			PsatCO2 <- 1.2264e7 * exp(-3167.8 / Ts); 
			Tsat <- -3167.8 / ln (pCO2 / 1.2264e7); // temperatura nasyceina - jaka powinna być wzgledem ciśnienia
			
			
			tCO2 <- 0.004 * ((pCO2)* Pa2bar / ga )^0.4551;
			tau <- tCO2 + tauDust;
		
			latentCO2 <- 0.0;

			insol_en <- ( (frozenCO2 > 1e-8)?(1 - albedoIce):(1 - albedoGlobal ) )* insol;  				// insolation ENERGY
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
				co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * (exp(exp_param * Tsat/Ts) - exp(exp_param))  ;
				
				
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
				co2_excess <- co2_excess_factor * totalCO2 * CO2freezingCoeff * (exp(exp_param * Ts/Tsat) - exp(exp_param)) ;
				//co2_excess <- CO2freezingCoeff * 10^10 * exp(-30000/(8.314462 * Ts))*44.0095/1000/ (marsArea / 4002);
				
				
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
			
			energy <- insol_en - rad_en + green_en + eddy_Ts + hadley_Ts  + latentCO2;  // latent heat
							   
			Ts <- prevTs			// previous step temp. 
				+  energy * delta_t / (cCO2       * pCO2 + cSoil     * soilDepth * soilDensity); // mass m-2 * conduction - without multiplication by area - its flux
//										  Wm-2 / J*kg-1*K-1 *  kg * m-2          + J*kg-1*K-1 * m         * kg * m-3
//										  Wm-2 1/ J*m-2*K-1                       + J*K-1*m-2
//										  Wm-2 /K-1 * m2 * J-1
//										  Ks-1 
		
		
			// check if co2_excess is too large
			float tmpPsatCO2 <- 1.2264e7 * exp(-3167.8 / Ts);
			float tmpPCO2 <- (atmCO2 + co2_excess * (pCO2greater?(1.0):(-1.0)) ) * exp(- height / 10000);
			
			if ( (tmpPCO2 > tmpPsatCO2 and pCO2greater) or (tmpPCO2 < tmpPsatCO2 and !pCO2greater)) {break; }
			else {co2_excess_factor <- co2_excess_factor * 0.5; }
		} 

		Tps <- Ts - Gamma_HT * height; // poprawka - zmiana znaku na "-" 2023-05-17
							
		heat_flux <- heat_flux - energy;	// obliczenie przepływu - różnicowo, pewnie trzeba dodać sol * delta_h, ale to później
		
			
	}
}
