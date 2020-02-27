#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>

using namespace std;

SPH_main domain;

void file_out(SPH_particle* part, int ite);

int main(void)
{
	//parse the input configuration file
	ifstream myfile("input.txt");
	string line;
	double t_total, t_print, delta_x;
	bool analyse_mode;
	vector<vector<double>> areas;
	vector<bool> particle_type;
	if (myfile.is_open())
	{
		getline(myfile, line);
		t_total = stod(line);
		getline(myfile, line);
		t_print = stod(line);
		getline(myfile, line);
		delta_x = stod(line);
		getline(myfile, line);
		analyse_mode = (line == "1");

		while (getline(myfile, line))
		{
			vector<double> this_area;
			size_t current, previous = 0;
			current = line.find(' ');
			while (current != string::npos) {
				this_area.push_back(stod(line.substr(previous, current - previous)));
				previous = current + 1;
				current = line.find(' ', previous);
			}
			areas.push_back(this_area);
			string type = line.substr(previous, line.length() - previous);
			particle_type.push_back(type == "1");
		}
	}
	else {
		cout << "Unable to open input configuration file" << endl;
		return -1;
	}

	domain.set_values(delta_x);									//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	for (int i = 0; i < areas.size(); i++)
		domain.place_points(areas[i][0], areas[i][1], areas[i][2], areas[i][3], particle_type[i]);

	domain.time_dynamic();
	domain.allocate_to_grid();

	stringstream name;
	name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << 0 << ".vtp";		
	write_file(name.str().c_str(), &domain.particle_list);

	double t_max = t_total;
	double t = 0;
	t += domain.dt;
	int count = 1;
	double dt_print = 0;
	
	while (t < t_max)
	{
		cout << "t = " << t << endl;
		// first half step
		#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
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
		#pragma omp parallel for schedule(static, 1) num_threads(4)
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);

		domain.reset_grid_count();

		#pragma omp parallel for schedule(static, 1) num_threads(4)
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.allocate_to_grid();								//update grid index of each particle
		
		// first full step
		#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
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
		#pragma omp parallel for schedule(static, 1) num_threads(4)
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.update_particle(&domain.particle_list[i]);
			domain.full_update(&domain.particle_list[i]);
		}
		
		#pragma omp parallel for schedule(static, 1) num_threads(4)
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.allocate_to_grid();								//update grid index of each particle
		domain.reset_grid_count();
		
		if (count % 10 == 0) {
			// cout << "Density field smoothed at iter = " << iter << endl;
			#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
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
			#pragma omp parallel for schedule(static, 1) num_threads(4)
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.update_rho(&domain.particle_list[i]);
		}
		domain.reset_grid_count();
		
		// get the max for minimum dynamic time step
		#pragma omp parallel for schedule(static, 1) num_threads(4)
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.get_new_max(&domain.particle_list[i]);

		// get new dynamic time step
		domain.time_dynamic();

		dt_print += domain.dt;
		if (dt_print >= t_print)
		{
			stringstream name;
			name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << (int)(t/t_print) << ".vtp";		
			write_file(name.str().c_str(), &domain.particle_list);
			dt_print = 0;
		}
		t += domain.dt;
		count++;
	}
	
	return 0;
}
