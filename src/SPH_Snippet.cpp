#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

SPH_main domain;

void file_out(SPH_particle* part, int ite);

int main(void)
{
	//parse the input configuration file
	ifstream myfile("input.txt");
	string line;
	double t_total, t_frame, delta_x;
	vector<vector<double>> areas;
	vector<bool> particle_type;
	if (myfile.is_open())
	{
		getline(myfile, line);
		t_total = stod(line);
		getline(myfile, line);
		t_frame = stod(line);
		getline(myfile, line);
		delta_x = stod(line);

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
		for (int i = 0; i < areas.size(); i++) {
			for (int j = 0; j < areas[i].size(); j++)
				cout << areas[i][j];
			cout << '\n' << particle_type[i] << endl;
		}
	}
	else cout << "Unable to open input configuration file" << endl;


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
	for (int iter = 1; iter < 10000; iter++) {

		// cout << "iter = " << iter << endl;
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
		domain.reset_grid_count();
		
		// get the max for minimum dynamic time step
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.get_new_max(&domain.particle_list[i]);

		// get new dynamic time step
		domain.time_dynamic();

		if (iter % 50 == 0)
		{
			stringstream name;
			name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << iter << ".vtp";		
			write_file(name.str().c_str(), &domain.particle_list);
		}
	}
	return 0;
}
