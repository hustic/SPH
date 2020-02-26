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
	domain.place_points(0.0, 10, 20.0, 10.52);				//top boundary
	domain.place_points(20.0, -0.52, 20.52, 10.52);				//right boundary
	domain.place_points(0.0, -0.52, 3.0, 5.0);
	domain.place_points(3.0, -0.52, 20.0, 2.0);
	domain.time_dynamic();
	domain.allocate_to_grid();

	stringstream name;
	name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << 0 << ".vtp";		
	write_file(name.str().c_str(), &domain.particle_list);

	double t_max = 30;
	double t = 0;
	t += domain.dt;
	int count = 1;
	double dt_print = 0;

	while (t < t_max)
	{
		// first half step
		for (int j = 0; j < domain.max_list[1]; j++)
		{
			for (int i = 0; i < domain.max_list[0]; i++)
			{
				for (int k = 0; k < domain.search_grid[i][j].size(); k++)
				{
					domain.store_initial(domain.search_grid[i][j][k]);
					domain.neighbour_iterate(domain.search_grid[i][j][k]);
				}
			}
		}
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);

		domain.reset_grid_count();

		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.allocate_to_grid();								//update grid index of each particle
		
		// first full step
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

		// last full step
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.update_particle(&domain.particle_list[i]);
			domain.full_update(&domain.particle_list[i]);

		}
		
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.allocate_to_grid();								//update grid index of each particle
		domain.reset_grid_count();
		
		if (count % 10 == 0) {
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
		if (count % 10 == 0) 
		{
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.update_rho(&domain.particle_list[i]);
		}
		domain.reset_grid_count();
		
		// get the max for minimum dynamic time step
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.get_new_max(&domain.particle_list[i]);

		// get new dynamic time step
		domain.time_dynamic();

		dt_print += domain.dt;
		if (dt_print >= 0.1)
		{
			stringstream name;
			name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << (int)(t/0.1) << ".vtp";		
			write_file(name.str().c_str(), &domain.particle_list);
			dt_print = 0;
		}
		t += domain.dt;
		count++;
	}
	return 0;
}
