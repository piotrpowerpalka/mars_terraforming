/**
* Name: testInsol
* Based on the internal empty template. 
* Author: drimn
* Tags: 
*/


model testInsol

/* Insert your model definition here */
global{
	float S0 <- 589;		// stała słoneczna
	float e <- 0.09341233;		// mimośród (spłaszczenie orbity)
	float nachylenieOsi <- 24.936;
	/**
	 * argument sin w GAMA jest podawany w stopniach
	 */
//	float Insol(float sollon, float lat) {
//		return Hobh / 24.0;
//	}
	init {
		loop sl from: 0 to: 360 step: 1 {
			loop lat from: -90 to: 90 step: 1  {
				create xxx with: [sollon::sl, lat::lat]{
					
					float fi <- lat;
					if (fi = 90.0)  { fi <- 89.99; } 
					if (fi = -90.0) { fi <- -89.99; }
					
					kat <- nachylenieOsi * sin(sollon);
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
					
					
					insol <- (24/#pi) * S0 * ((1 + e * cos((sollon - 248) mod 360))^2) / ((1 - e^2)^2) 
								* max(0.0, sin(fi) * sin(nachylenieOsi) * sin(sollon) 
								* OMEGA * 2 * #pi / 360.0 + cos(fi) * cos(kat) * sin(OMEGA)
								) / 24;
		
					write "" + sl + ";" + lat + ";" + insol;
					
				}
			}
		}
	}

	species xxx {
		float sollon;
		float lat;
		float insol;
		
		float kat;
		float zachSlonca;
		float OMEGA;
		
		aspect base {
			
			draw square(1.0) at: {sollon/2, lat} color: hsb(insol/1000, 1, 1); 
		}		
	}
}
experiment Matrices {
	output {
		display mars_plain_co2 type: opengl ambient_light: 100 background: #white  {
			species xxx aspect: base;
		}
	}
}

