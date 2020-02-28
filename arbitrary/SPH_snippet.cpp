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
	ifstream myfile("arbitrary.txt");
	string line;
	double t_total, t_print, delta_x;
	if (myfile.is_open())
	{
		getline(myfile, line);
		t_total = stod(line);
		getline(myfile, line);
		t_print = stod(line);
		getline(myfile, line);
		delta_x = stod(line);
		//domain.set_values(delta_x);
		//domain.initialise_grid();
		getline(myfile, line);
		domain.vcount = stoi(line);
		domain.vertices.resize(domain.vcount);
		domain.boundaries.resize(domain.vcount);
		for (int i = 0; i < domain.vcount; i++)
		{
			domain.vertices[i].resize(2);
			domain.boundaries[i].resize(4);
			getline(myfile, line);
			size_t current, previous = 0;
			current = line.find(' ');
			domain.vertices[i][0] = stod(line.substr(previous, current - previous));
			domain.vertices[i][1] = stod(line.substr(current, line.length() - current));
		}
		for (int i = 0; i < domain.vcount; i++)
		{
			domain.boundaries[i][0] = domain.vertices[i][0];
			domain.boundaries[i][1] = domain.vertices[i][1];
			domain.boundaries[i][2] = domain.vertices[(i + 1) % (domain.vcount)][0];
			domain.boundaries[i][3] = domain.vertices[(i + 1) % (domain.vcount)][1];
		}
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
			this_area.push_back(stod(line.substr(previous, line.length() - previous)));
			domain.areas.push_back(this_area);
		}
	}
	else {
		cout << "Unable to open input configuration file" << endl;
		return -1;
	}

	
	domain.set_values(delta_x);										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	for (int i = 0; i < domain.vcount; i++)
	{
		domain.place_points_new(domain.boundaries[i][0], domain.boundaries[i][1], domain.boundaries[i][2], domain.boundaries[i][3], 1);
	}
	for (int i = 0; i < domain.areas.size(); i++)
	{

		domain.place_points_new(domain.areas[i][0], domain.areas[i][1], domain.areas[i][2], domain.areas[i][3], 0);
	}

	//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	//domain.place_points_new(0.0, 10.3, 0.0, -0.31, 1);
	//domain.place_points_new(19.9, 10.0, 0.05, 10.0, 1);
	//domain.place_points_new(10.0, 0.0, 20.0, 10.0, 1);
	//domain.place_points_new(0.1, 0.0, 10.0, 0.0, 1);
	//domain.place_points_new(0.1, 0.1, 3.05, 5.05, 0);
	//domain.place_points_new(3.1, 0.1, 19.95, 2.05, 0);

	//domain.place_points_new(0.0, 10.6, 0.0, -0.62, 1);
	//domain.place_points_new(19.8, 10.0, 0.1, 10.0, 1);
	//domain.place_points_new(10.0, 0.0, 20.0, 10.0, 1);
	//domain.place_points_new(0.2, 0.0, 10.0, 0.0, 1);
	//domain.place_points_new(0.2, 0.2, 3.1, 5.1, 0);
	//domain.place_points_new(3.2, 0.2, 19.9, 2.1, 0);
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
		domain.time_dynamic();
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.reset_grid_count();
		domain.allocate_to_grid();								//update grid index of each particle

		dt_print += domain.dt;
		if (dt_print >= t_print)
		{
			cout << "t = " << t << endl;
			stringstream name;
			name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << (int)(t / t_print) << ".vtp";
			write_file(name.str().c_str(), &domain.particle_list);
			dt_print = 0;
		}
		t += domain.dt;
		count++;
	}

	return 0;

}