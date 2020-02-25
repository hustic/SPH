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


	//write_file(name.str().c_str(), &domain.particle_list);{

	
	for (int iter = 1; iter < 30; iter++) {
		cout << "iter = " << iter << endl;
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			if (domain.particle_list[i].x[0] < domain.min_x[0] || domain.particle_list[i].x[1] < domain.min_x[1])
				cout << "out of min range" << endl;
			else if (domain.particle_list[i].x[0] > domain.max_x[0] || domain.particle_list[i].x[1] > domain.max_x[1])
				cout << "out of max range" << endl;

				domain.neighbour_iterate(&domain.particle_list[i]);
		}
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);
		if (iter % 10 == 0) {
			cout << "Density field smoothed at iter = " << iter << endl;
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.density_field_smoothing(&domain.particle_list[i]);
		}
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.allocate_to_grid();								//update grid index of each particle
		
		stringstream name;
        name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << iter << ".vtp";
		

		//write_file(name.str().c_str(), &domain.particle_list);
		
	}
	

	
	
	return 0;
}
