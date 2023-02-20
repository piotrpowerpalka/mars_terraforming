/**
* Name: basic
* Based on the internal empty template. 
* Author: pawel
* Tags: 
*/


model basic

/* Insert your model definition here */

global {
	
	init {

    	create cell from: csv_file("../includes/outpkty_z_neigh_4002_fit.csv",";",true) with:
    	[
    		id_cell::int(read("id")),
    		x::int(read("x")),
	    	y::int(read("y")),
    		z::int(read("z")),     	
    		longitude::float(read("longitude")),
   		 	latitude::float(read("lattitude")),
    		spec::int(read("spec")),
    		
    		n[0]::int(read("n1")),
    		n[1]::int(read("n2")),
    		n[2]::int(read("n3")),
    		n[3]::int(read("n4")),
    		n[4]::int(read("n5")),
    		n[5]::int(read("n6")),
    		
    		a[0]::int(read("a1")),
    		a[1]::int(read("a2")),
    		a[2]::int(read("a3")),
    		a[3]::int(read("a4")),
    		a[4]::int(read("a5")),
    		a[5]::int(read("a6")),
    		
    		temp::float(read("temperture")),
    		zonal_wind::float(read("zonal_wind")),
    		meridional_wind::float(read("meridional_wind"))
    		
    	];
    	
    }
    reflex initialize when: cycle = 0 {
						
			if (spec = 1) {
				remove index: 6 from: n;
				remove index: 6 from: a;
			}
			neighbors <- cell where (each.id_cell in n);
			
			int ntmp;
		int atmp; 
		
		loop idx from: 0 to: n.length { 
			loop jdx from: 0 to: n.length-1 {
				if (a[idx] > a[jdx]) {
					atmp  <- a[idx];
					a[idx] <- a[jdx];
					a[jdx] <- atmp;
					
					ntmp  <- n[idx];
					n[idx] <- n[jdx];
					n[jdx] <- ntmp;
				}
			}
		}
			
			if (spec != 1) {
    		
    		sh <- polygon([{(neighbors[0].longitude + longitude)/2, (neighbors[0].latitude + latitude)/2}, 
				{(neighbors[1].longitude + longitude)/2, (neighbors[1].latitude + latitude)/2},
				{(neighbors[2].longitude + longitude)/2, (neighbors[2].latitude + latitude)/2},
				{(neighbors[3].longitude + longitude)/2, (neighbors[3].latitude + latitude)/2},
				{(neighbors[4].longitude + longitude)/2, (neighbors[4].latitude + latitude)/2},
				{(neighbors[5].longitude + longitude)/2, (neighbors[5].latitude + latitude)/2}
			]);
    		
    	} else {
    		sh <- polygon([{(neighbors[0].longitude + longitude)/2, (neighbors[0].latitude + latitude)/2}, 
				{(neighbors[1].longitude + longitude)/2, (neighbors[1].latitude + latitude)/2},
				{(neighbors[2].longitude + longitude)/2, (neighbors[2].latitude + latitude)/2},
				{(neighbors[3].longitude + longitude)/2, (neighbors[3].latitude + latitude)/2},
				{(neighbors[4].longitude + longitude)/2, (neighbors[4].latitude + latitude)/2}
			]);
    	} 
			
			 
	    }
    
}



experiment basic type: gui {
	output{
		display mars type: opengl ambient_light: 100
		{
			//species cell aspect: base;
		}
		
	}
}
