﻿#include "SPH_2D.h"
#include <omp.h>


SPH_main *SPH_particle::main_data;


SPH_main::SPH_main(): grid_count(), v_max(v_max), a_max(a_max), rho_max(rho_max)
{
	SPH_particle::main_data = this;
}


//function to calculate which grid each cell belongs to
void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
	{
		if (x[i] < main_data->min_x[i])
		{
			x[i] = main_data->min_x[i]; 
		}
		else if (x[i] > main_data->max_x[i])
		{
			x[i] = main_data->max_x[i];
		}
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
	}
}



double SPH_main::cubic_spline(double r[2])
{
	double q = sqrt(r[0] * r[0] + r[1] * r[1]) / h;
	if (q >= 0 && q <= 1)
	{
		return 10 * (1 - 1.5 * q * q + 0.75 * pow(q, 3)) / (7 * M_PI * h * h);
	}
	else if (q > 1 and q <= 2)
	{
		return 10  * 0.25 * pow((2 - q), 3) / (7 * M_PI * h * h);
	}
	else { return 0; }
}

double SPH_main::cubic_spline_first_derivative(double r[2])
{
	double q = sqrt(r[0] * r[0] + r[1] * r[1]) / h;
	if (q >= 0 && q <= 1)
	{
		return 10 * (-3 * q + 2.25 * q * q) / (7 * M_PI * pow(h, 3));
	}
	else if (q > 1 and q <= 2)
	{
		return -10 * 0.75 * (2 - q) * (2 - q) / (7 * M_PI * pow(h, 3));
	}
	else { return 0; }
}

void SPH_main::update_gradients(double r[2], SPH_particle* part, SPH_particle* other_part)		//updates acceleration and rate of change of density of a particle
{	
	
	double vij[2], eij[2];
	double dwdr = cubic_spline_first_derivative(r);
	for (int n = 0; n < 2; n++)
	{
		vij[n] = part->v[n] - other_part->v[n];
		eij[n] = r[n] / sqrt(r[0] * r[0] + r[1] * r[1]);
		#pragma omp atomic
		part->a[n] += -mass * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dwdr * eij[n] + mu * mass * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dwdr * vij[n] / sqrt(r[0] * r[0] + r[1] * r[1]);

	}
	if (part->vij_half[0] < sqrt(vij[0] * vij[0] + vij[1] * vij[1]) && (!part->is_boundary || !other_part->is_boundary))
	{
		part->vij_half[0] = sqrt(vij[0] * vij[0] + vij[1] * vij[1]);
	}

	#pragma omp atomic
	part->D += mass * dwdr * (vij[0] * eij[0] + vij[1] * eij[1]);
}


void SPH_main::density_field_smoothing(SPH_particle* part)		//performs the density field smoothing for a particle
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn1[2];			//vector from 1st to 2nd particle
	double dn2[2];			//vector from 1st to 2nd particle

	int cnt;
	// #pragma parallel for
	for (int j = part->list_num[1]; j <= part->list_num[1] + 1; j++)
		if (j >= 0 && j < max_list[1])
			for (int i = part->list_num[0] - j + part->list_num[1]; i <= part->list_num[0] + 1; i++)
				if (i >= 0 && i < max_list[0])
				{
					// if not in the same grid
					if (j != part->list_num[1] || i != part->list_num[0])
					{
						for (cnt = 0; cnt < search_grid[i][j].size(); cnt++)
						{
							other_part = search_grid[i][j][cnt];

							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								dn1[n] = part->x[n] - other_part->x[n];
								dn2[n] = other_part->x[n] - part->x[n];
							}

							dist = sqrt(dn1[0] * dn1[0] + dn1[1] * dn1[1]);
							if (dist < 2. * h)		//only particle within 2h
							{
								#pragma omp atomic
								part->numerator += cubic_spline(dn1);
								#pragma omp atomic
								part->denominator += cubic_spline(dn1) / other_part->rho;

								#pragma omp atomic
								other_part->numerator += cubic_spline(dn2);
								#pragma omp atomic
								other_part->denominator += cubic_spline(dn2) / part->rho;

								if (analysisMode)
								{
									part->rho_spacing += cubic_spline(dn1) * mass / rho0;
									other_part->rho_spacing += cubic_spline(dn2) * mass / rho0;
								}
							}
						}
					}
					// if in the same grid
					else
					{
						for (cnt = grid_count[i][j]; cnt < search_grid[i][j].size(); cnt++)
						{
							other_part = search_grid[i][j][cnt];

							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								dn1[n] = part->x[n] - other_part->x[n];
								dn2[n] = other_part->x[n] - part->x[n];
							}
							dist = sqrt(dn1[0] * dn1[0] + dn1[1] * dn1[1]);
							if (dist < 2. * h)					//only particle within 2h
							{
								#pragma omp atomic
								part->numerator += cubic_spline(dn1);
								#pragma omp atomic
								part->denominator += cubic_spline(dn1) / other_part->rho;

								#pragma omp atomic
								other_part->numerator += cubic_spline(dn2);
								#pragma omp atomic
								other_part->denominator += cubic_spline(dn2) / part->rho;

								if (analysisMode)
								{
									part->rho_spacing += cubic_spline(dn1) * mass / rho0;
									other_part->rho_spacing += cubic_spline(dn2) * mass / rho0;
								}
							}
						}
						++grid_count[i][j];
					}
				}
	if (part->numerator != 0) 
	{
		part->rho2 = part->numerator / part->denominator;
	}
}

// set initial parameters
void SPH_main::set_values(double delta_x)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = delta_x;
	c0 = 20;

	mu = 0.001;
	g[0] = 0.0;
	g[1] = -9.81;
	rho0 = 1000;
	mass = rho0 * dx * dx;
	h_fac = 1.3;

	cfl = 0.1;

	a_max = -g[1];
	v_max = 0;
	rho_max = rho0;

	h = dx*h_fac;
	// dt = 0.5 * 0.1 * h / c0;

	gamma = 7.0;
	B = c0 * c0 * rho0 / gamma;

}

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 3.0 * dx;
		max_x[i] += 3.0 * dx; //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	grid_count.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
	{
		search_grid[i].resize(max_list[1]);
		grid_count[i].resize(max_list[1]);
	}
}

// place particles into the grid
void SPH_main::place_points(double min0, double min1, double max0, double max1, bool type)
{
	double x[2] = { min0, min1 };
	SPH_particle particle;

	while (x[0] <= max0)
	{
		x[1] = min1;
		while (x[1] <= max1)
		{
			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.a[i] = 0.0 + g[i];
				particle.v[i] = 0.0;
			}
			particle.D = 0.0;
			particle.rho = rho0;
			particle.P = 0.0;

			particle.is_boundary = type;
			
			particle.calc_index();
			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}
}


void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
	#pragma omp parallel for schedule(static, 1) num_threads(4)
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
	a_max = 0;
	v_max = 0;
	rho_max = 0;
}



void SPH_main::neighbour_iterate(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn1[2];			//vector from 1st to 2nd particle
	double dn2[2];
    int cnt;


	for (int j = part->list_num[1]; j <= part->list_num[1] + 1; j++)
		if (j>= 0 && j < max_list[1])
			for (int i = part->list_num[0] - j + part->list_num[1]; i <= part->list_num[0] + 1; i++)
				if (i >= 0 && i < max_list[0])
				{
					// if not in the same grid
					if (j != part->list_num[1] || i != part->list_num[0])
					{
						for (cnt = 0; cnt < search_grid[i][j].size(); cnt++)
						{
							other_part = search_grid[i][j][cnt];

							if (part != other_part)		//stops particle interacting with itself
							{
								//Calculates the distance between potential neighbours
								for (int n = 0; n < 2; n++)
								{
									dn1[n] = part->x[n] - other_part->x[n];
									dn2[n] = other_part->x[n] - part->x[n];
								}

								dist = sqrt(dn1[0] * dn1[0] + dn1[1] * dn1[1]);
								if (dist < 2. * h && dist != 0)					//only particle within 2h
								{
									update_gradients(dn1, part, other_part);
									update_gradients(dn2, other_part, part);
								}
							}
						}
					}
					// if in the same grid
					else
					{
						for (cnt = grid_count[i][j]; cnt < search_grid[i][j].size(); cnt++)
						{
							other_part = search_grid[i][j][cnt];

							if (part != other_part)		//stops particle interacting with itself
							{
								//Calculates the distance between potential neighbours
								for (int n = 0; n < 2; n++)
								{
									dn1[n] = part->x[n] - other_part->x[n];
									dn2[n] = other_part->x[n] - part->x[n];
								}

								dist = sqrt(dn1[0] * dn1[0] + dn1[1] * dn1[1]);
								if (dist < 2. * h && dist != 0)					//only particle within 2h
								{
									update_gradients(dn1, part, other_part);
									update_gradients(dn2, other_part, part);
								}
							}
						}
						++grid_count[i][j];
					}
				}

}

// update the particle with half step
void SPH_main::update_particle(SPH_particle* part) 
{
	if (!part->is_boundary)
	{
		for (int k = 0; k < 2; k++)
		{
			if (part->x[k] < min_x[k] + 3.5 * dx)
			{
				double dist = abs(part->x[k] - (min_x[k] + 3.0 * dx));
				#pragma omp atomic
				part->a[k] += abs(repulsion(part, dist));
			}
			else if(part->x[k] > max_x[k] - 3.5 * dx){

				double dist = abs(part->x[k] - (min_x[k] + 3.0 * dx));
				#pragma omp atomic
				part->a[k] -= abs(repulsion(part, dist));
			}
            part->x[k] = part->x_half[k] + 0.5 * dt * part->v[k];
            part->v[k] = part->v_half[k] + 0.5 * dt * part->a[k];
		}
		for (int n = 0; n < 2; n++)
		{
			part->a[n] = 0.0 + g[n];
		}
	}

	part->rho = part->rho_half + 0.5 * part->D * dt;
	part->P = B * (pow((part->rho / rho0), gamma) - 1);

	part->D = 0.0;
}


void SPH_main::reset_grid_count()
{
	for (auto& temp : grid_count) fill(temp.begin(), temp.end(), 0);

}

// update rho for density smoothing and reset some parameters
void SPH_main::update_rho(SPH_particle* part)
{
	if (part->rho2 != 0) 
	{
		part->rho = part->rho2;
	}
	part->rho2 = 0;
	part->numerator = 0.0;
	part->denominator = 0.0;
}

// store the initial velocity density and position for second order scheme
void SPH_main::store_initial(SPH_particle* part)
{
	part->rho_half = part->rho;
	part->v_half[0] = part->v[0];
	part->v_half[1] = part->v[1];
	part->x_half[0] = part->x[0];
	part->x_half[1] = part->x[1];

}

// the final update of particle for each iteration using the second order method
void SPH_main::full_update(SPH_particle* part)
{
	if (!part->is_boundary)
	{
		for (int k = 0; k < 2; k++)
		{
            part->x[k] = 2 * part->x[k] - part->x_half[k];
            part->v[k] = 2 * part->v[k] - part->v_half[k];
			
			part->a_half[k] = part->a[k];
		}
		for (int n = 0; n < 2; n++)
		{
			part->a[n] = 0.0 + g[n];
		}
	}

	part->rho = 2 * part->rho - part->rho_half;
	part->P = B * (pow((part->rho / rho0), gamma) - 1);

	part->D = 0.0;
}

// Find the timestep for next iteration by finding the minimum of three other timesteps
// calculated using the maximum velocity, density and acceleration.
void SPH_main::time_dynamic()
{
	dt_cfl = (h / v_max);
	dt_a = sqrt(h / a_max);
	dt_f = h / (c0 * sqrt(pow((rho_max / rho0), gamma - 1)));

	if (dt_cfl <= dt_a && dt_cfl <= dt_f && dt_cfl != 0)
	{
		dt = cfl * dt_cfl;
	}
	else if (dt_a <= dt_cfl && dt_a <= dt_f && dt_a != 0)
	{
		dt = cfl * dt_a;
	}
	else if (dt_a <= dt_f && dt_a != 0 && dt_cfl != 0)
	{
		dt = cfl * dt_a;
	}
	else if (dt_f != 0)
	{
		dt = cfl * dt_f;
	}
	else
	{
		dt = 0.5 * 0.1 * h / c0;
	}
	dt_f = 0;
	dt_a = 0;
	dt_cfl = 0;
}

// get the new global maximum velocity acceleration and density
void SPH_main::get_new_max(SPH_particle* part)
{
	double temp_vv = 0;
	double temp_aa = 0;

	if (!part->is_boundary) 
	{
		temp_vv = part->vij_half[0];
		temp_aa = sqrt(part->a_half[0] * part->a_half[0] + part->a_half[1] * part->a_half[1]);
		if (v_max < temp_vv) v_max = temp_vv;
		if (a_max < temp_aa) a_max = temp_aa;
		if (rho_max < part->rho) rho_max = part->rho;
	}
}

// calculate the repulsive force to the boundary
double SPH_main::repulsion(SPH_particle *part, double &dist)
{
	double fd = pow((0.5 * dx / dist), 6) - 1.0;
	double temp_a = fd * 9.81 * 1.0 / dx;

    return temp_a;
}

// update function for the first order method
void SPH_main::update_particle_FE(SPH_particle* part) 
{
	if (!part->is_boundary)
	{
		for (int k = 0; k < 2; k++)
		{
			if (part->x[k] < min_x[k] + 3.5 * dx)
			{
				double dist = abs(part->x[k] - (min_x[k] + 3.0 * dx));
				#pragma omp atomic
				part->a[k] += abs(repulsion(part, dist));
			}
			else if(part->x[k] > max_x[k] - 3.5 * dx){

				double dist = abs(part->x[k] - (min_x[k] + 3.0 * dx));
				#pragma omp atomic
				part->a[k] -= abs(repulsion(part, dist));
			}
            part->x[k] = part->x[k] + dt * part->v[k];
            part->v[k] = part->v[k] + dt * part->a[k];
		}
		for (int n = 0; n < 2; n++)
		{
			part->a[n] = 0.0 + g[n];
		}
	}

	part->rho = part->rho + part->D * dt;
	part->P = B * (pow((part->rho / rho0), gamma) - 1);

	part->D = 0.0;
}