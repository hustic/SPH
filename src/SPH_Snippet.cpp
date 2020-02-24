#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

SPH_main domain;

int main(void)
{
	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x,domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

	domain.allocate_to_grid();									//needs to be called for each time step

	stringstream name;
	name << "initial_configuration.vtp";


	write_file(name.str().c_str(), &domain.particle_list);


	for (int iter = 0; iter < 10; iter++) {
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.neighbour_iterate(&domain.particle_list[i]);
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);
		if (iter % 10 == 0) {
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.density_field_smoothing(&domain.particle_list[i]);
		}
		
		
		stringstream name;
        name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << iter << ".vtp";
		

		write_file(name.str().c_str(), &domain.particle_list);

	}
	

	
	
	return 0;
}
