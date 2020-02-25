#include "SPH_2D.h"
#include "file_writer.h"
#include "gnuplot.h"
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
	domain.place_points(0.0, 5.0, 3.52, 10.52);
	domain.place_points(3.0, 7.0, 20.0, 10.52);

	domain.allocate_to_grid();

	fstream fp;
	stringstream name;
	name << "ite" << "_" << 0 << ".dat";
	fp.open(name.str().c_str(), ios_base::out);
	for (int i = 0; i < domain.particle_list.size(); i++)
	{
		fp << domain.particle_list[i].x[0] << " " << domain.particle_list[i].x[1] << "\n";
	}
	fp.close();
	// name << "initial_configuration.vtp";


	// write_file(name.str().c_str(), &domain.particle_list);

	for (int iter = 1; iter < 10000; iter++) {

		// cout << "iter = " << iter << endl;
		for (int i = 0; i < domain.particle_list.size(); i++)
		{
			
			// if (domain.particle_list[i].x[0] < domain.min_x[0] || domain.particle_list[i].x[1] < domain.min_x[1])
				// cout << "out of min range" << endl;
			// else if (domain.particle_list[i].x[0] > domain.max_x[0] || domain.particle_list[i].x[1] > domain.max_x[1])
				// cout << "out of max range" << endl;

				domain.neighbour_iterate(&domain.particle_list[i]);
		}
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_particle(&domain.particle_list[i]);

		domain.reset_grid_count();
		if (iter % 10 == 0) {
			// cout << "Density field smoothed at iter = " << iter << endl;
			for (int i = 0; i < domain.particle_list.size(); i++)
				domain.density_field_smoothing(&domain.particle_list[i]);
		}
		if (iter % 10 == 0) 
		{
		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.update_rho(&domain.particle_list[i]);
		}

		// cout << domain.particle_list[2000].x[0] << " " << domain.particle_list[2000].x[1] << endl;

		for (int i = 0; i < domain.particle_list.size(); i++)
			domain.particle_list[i].calc_index();

		domain.reset_grid_count();
		domain.allocate_to_grid();								//update grid index of each particle
		stringstream name;
        // name << "output" << "_" << setfill('0') << setw(int(to_string(100).length())) << iter << ".vtp";		
		
		if (iter % 100 == 0)
		{
			name << "ite" << "_" << iter << ".dat";
			fp.open(name.str().c_str(), ios_base::out);
			// write_file(name.str().c_str(), &domain.particle_list);
			for (int i = 0; i < domain.particle_list.size(); i++)
			{

				fp << domain.particle_list[i].P << " " << domain.particle_list[i].rho << "\n";
			}
			fp.close();
		}
	}
	return 0;
}
