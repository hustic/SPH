#include "SPH_2D.h"
#include "file_writer.h"
#include <string>

SPH_main domain;

int main(void)
{
	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x,domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

	domain.allocate_to_grid();									//needs to be called for each time step


	for (int iter = 1; iter < 30; iter++) {
		domain.neighbour_iterate(&domain.particle_list[100]);		//finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle
		if (iter % 10 == 0) {
			domain.density_field_smoothing(&domain.particle_list[100]);
		}
		
		char name [7];
		string str = to_string(iter);
		for (int i = 0; i < str.length(); i++)
			str[2 - i] = name[i];
		for (int i = 0; i < 3 - str.length(); i++)
			str[i] = '0';

		write_file(name, &domain.particle_list);

	}
	

	
	
	return 0;
}
