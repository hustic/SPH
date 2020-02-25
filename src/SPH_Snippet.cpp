#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

SPH_main domain;

void file_out(SPH_particle* part, int ite);

int main(void)
{
	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	domain.place_points(-0.52, -0.52, 0.0, 10.52);				//left boundary
	domain.place_points(0.0, -0.52, 20.0, 0.0);				//top boundary
	domain.place_points(20.0, -0.52, 20.52, 10.52);				//right boundary
	domain.place_points(0.0, 5.0, 3.0, 10.52);
	domain.place_points(3.0, 8.0, 20.0, 10.52);
	domain.allocate_to_grid();

	stringstream name;
	name << "initial_configuration.vtp";
	write_file(name.str().c_str(), &domain.particle_list);

	for (int iter = 1; iter < 5000; iter++) {

		// cout << "iter = " << iter << endl;
		for (int j = 0; j < domain.max_list[1]; j++)
		{
			for (int i = 0; i < domain.max_list[0]; i++)
			{
				for (int k = 0; k < domain.search_grid[i][j].size(); k++)
				{
					domain.neighbour_iterate(domain.search_grid[i][j][k]);
				}
			}
		}

		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);

		domain.reset_grid_count();
		if (iter % 10 == 0) {
			// cout << "Density field smoothed at iter = " << iter << endl;
			for (int j = 0; j < domain.max_list[1]; j++)
			{
				for (int i = 0; i < domain.max_list[0]; i++)
				{
					for (int k = 0; k < domain.search_grid[i][j].size(); k++)
					{
						domain.density_field_smoothing(domain.search_grid[i][j][k]);
					}
				}
			}
		}
		if (iter % 10 == 0) 
		{
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_rho(&domain.particle_list[i]);
		}

		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.reset_grid_count();
		domain.allocate_to_grid();								//update grid index of each particle
		
		if (iter % 50 == 0)
		{
			stringstream name;
			name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << iter << ".vtp";		
			write_file(name.str().c_str(), &domain.particle_list);
		}
	}
	return 0;
}
