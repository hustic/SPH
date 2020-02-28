
#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class SPH_main;

class SPH_particle
{
public:
	double x[2], v[2];				//position and velocity
	double rho, P;					//density and pressure
	double a[2];					//acceleration
	double D;						//rate of change of density
	double rho2 = 0;
	double vij_half[2];

	double x_half[2], v_half[2];		//position and velocity
	double rho_half;					//density and pressure
	double a_half[2];					//acceleration
	double D_half;						//rate of change of density

	double numerator = 0;
	double denominator = 0;

	static SPH_main* main_data;		//link to SPH_main class so that it can be used in calc_index

	int list_num[2];				//index in neighbour finding array
	bool is_boundary;				//if it is a boundary particle
	void calc_index(void);

	double rho_spacing = 0;				//spacing density of particle for analysis
};


class SPH_main
{
public:
	SPH_main();
	double cubic_spline(double r[2]);
	double cubic_spline_first_derivative(double r[2]);
	void update_gradients(double r[2], SPH_particle* part, SPH_particle* other_part);
	void density_field_smoothing(SPH_particle* part);

	void set_values(double delta_x);
	void initialise_grid(void);

	void place_points(double min0, double min1, double max0, double max1, bool type);

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void neighbour_iterate(SPH_particle* part);
	void update_particle(SPH_particle* part);
	void reset_grid_count();
	void update_rho(SPH_particle* part);
	void store_initial(SPH_particle* part);
	void time_dynamic();
	void full_update(SPH_particle* part);
	void get_new_max(SPH_particle* part);
	double repulsion(SPH_particle* part, double &dist);
	void update_particle_FE(SPH_particle* part);
    
	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing
	double c0;								//speed of sound

	double dt;								//time step

	double g[2];							//gravity constant
	double mu;								//viscosity
	double rho0;							//initial density
	double B;
	double gamma;
	double mass;

	// for dynamic time stepping
	double v_max;
	double a_max;
	double rho_max;
	double dt_cfl;
	double dt_f;
	double dt_a;
	double cfl;

	double min_x[2], max_x[2];				//dimensions of simulation region

	vector<vector<int>> grid_count;

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell

	bool analysisMode;
};
