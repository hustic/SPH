#include "SPH_2D.h"
// #include <omp.h>

SPH_main *SPH_particle::main_data;


SPH_main::SPH_main(): grid_count(), v_max(v_max), a_max(a_max), rho_max(rho_max)
{
	SPH_particle::main_data = this;
}


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
	if (is_boundary) {
		main_data->near_boundary.insert({ list_num[0], list_num[1] });
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
		part->a[n] += -mass * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dwdr * eij[n] + mu * mass * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dwdr * vij[n] / sqrt(r[0] * r[0] + r[1] * r[1]);
	}
	/*
	if (part->vij_half[0] < sqrt(vij[0] * vij[0] + vij[1] * vij[1]) && !part->is_boundary && !other_part->is_boundary)
	{
		part->vij_half[0] = sqrt(vij[0] * vij[0] + vij[1] * vij[1]);
	}
	*/
	
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
								part->numerator += cubic_spline(dn1);
								part->denominator += cubic_spline(dn1) / other_part->rho;

								other_part->numerator += cubic_spline(dn2);
								other_part->denominator += cubic_spline(dn2) / part->rho;
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
								part->numerator += cubic_spline(dn1);
								part->denominator += cubic_spline(dn1) / other_part->rho;

								other_part->numerator += cubic_spline(dn2);
								other_part->denominator += cubic_spline(dn2) / part->rho;

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

void SPH_main::set_values(double delta_x)
{
	vector<double> values;
	values.resize(vcount);
	for (int i = 0; i < vcount; i++)
		values[i] = vertices[i][0];
	min_x[0] = *min_element(values.begin(), values.end()) - 1.0;
	max_x[0] = *max_element(values.begin(), values.end()) + 1.0;
	
	for (int i = 0; i < vcount; i++)
		values[i] = vertices[i][1];
	min_x[1] = *min_element(values.begin(), values.end()) - 1.0;
	max_x[1] = *max_element(values.begin(), values.end()) + 1.0;
	//min_x[0] = -1.0;
	//min_x[1] = -1.0;

	//max_x[0] = 21.0;
	//max_x[1] = 11.0;

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
	//boundaries.resize(4);
	//for (int i = 0; i < boundaries.size(); i++)
	//	boundaries[i].resize(4);
	//boundaries[0] = { 0.0, 10.0, 0.0, 0.0 };
	//boundaries[1] = { 0.0, 0.0, 9.98, 0.0 };
	//boundaries[2] = { 9.98, 0.0, 19.98, 10.0};
	////boundaries[3] = { 20.0, 5.1, 20.0, 10.0 };
	//boundaries[3] = { 19.98, 10.0, 0.0, 10.0 };

}

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		//min_x[i] -= 3.0 * dx;
		//max_x[i] += 3.0 * dx; //add buffer for virtual wall particles

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
			// cout << "Particle has been placed at (" << particle.x[0] << " ," << particle.x[1] << ") with type = " << particle.is_boundary<<endl;

			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
	}
}

int SPH_main::is_inside(double x, double y)
{
	
	
	for (int i = 0; i < boundaries.size(); i++)
	{
		double cross_product = (x - boundaries[i][0]) * (y - boundaries[i][3]) - (x - boundaries[i][2]) * (y - boundaries[i][1]);
		if (cross_product <= 0) {
			return i + 1;
		}
	}
	return 0;
}

void SPH_main::place_points_new(double x0, double y0, double x1, double y1, bool type)
{
	double x[2] = { x0, y0 };
	SPH_particle particle;
	if (type) {
		double slope;
		double increment[2];
		double len = sqrt((y1 - y0) * (y1 - y0) + (x1 - x0) * (x1 - x0));
		int cnt = (int)(len / dx) + 1;
		if (x0 != x1) {
			slope = (y1 - y0) / (x1 - x0);
			increment[0] = (x1 - x0) * dx / len;
			increment[1] = increment[0] * slope;
		}
		else {
			increment[0] = 0.0;
			if (y1 > y0) {
				increment[1] = dx;
			} else { increment[1] = -dx; }
		}
		for (int layer = 0; layer <= 3; layer++) {
			x[0] = x0+ increment[1] * layer;
			x[1] = y0 - increment[0] * layer;
			for(int c = 0; c < cnt; c++)
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
				particle.is_boundary = true;
				particle.calc_index();

				particle_list.push_back(particle);

				x[0] += increment[0];
				x[1] += increment[1];
			}
		}
	}else {
		x[0] = min(x0, x1);
		while (x[0] <= max(x0, x1))
		{
			x[1] = min(y0, y1);
			while (x[1] <= max(y0, y1))
			{
				if (is_inside(x[0], x[1])==0) 
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
					particle.is_boundary = false;
					particle.calc_index();
					// cout << "Particle has been placed at (" << particle.x[0] << " ," << particle.x[1] << ") with type = " << particle.is_boundary<<endl;

					particle_list.push_back(particle);
				}
				x[1] += dx;
			}
			x[0] += dx;
		}
	}

	
}

void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
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

void SPH_main::update_particle(SPH_particle* part) 
{
	if (!part->is_boundary  )
	{
		//double dist;
		if (near_boundary.find({ part->list_num[0],part->list_num[1] }) != near_boundary.end() && is_inside(part->x[0] + dt * part->v[0], part->x[1] + dt * part->v[1])!=0) {
			int outnum = is_inside(part->x[0] + dt * part->v[0], part->x[1] + dt * part->v[1]);
			//if (outnum != 0) {
				if (boundaries[outnum-1][0] == boundaries[outnum-1][2]) //out of vertical boundary
				{
					part->v[0] = -part->v[0];
				}
				else if (boundaries[outnum-1][1] == boundaries[outnum-1][3]) //out of horizontal boundary
				{
					part->v[1] = -part->v[1];
				}
				else // out of tilted boundary
				{
					double b[2], diff;
					b[0] = boundaries[outnum - 1][2] - boundaries[outnum - 1][0];
					b[1] = boundaries[outnum - 1][3] - boundaries[outnum - 1][1];
					diff = (part->v[0] * b[0] + part->v[1] * b[0]) / sqrt(b[0] * b[0] + b[1] * b[1]);
					part->v[0] = part->v[0] - 2. * (part->v[0] - diff * b[0] / sqrt(b[0] * b[0] + b[1] * b[1]));
					part->v[1] = part->v[1] - 2. * (part->v[1] - diff * b[1] / sqrt(b[0] * b[0] + b[1] * b[1]));
					
					/*double k = part->v[0];
					part->v[0] = part->v[1];
					part->v[1] = k;*/
				}
			//}
			
		}
		else {
			for (int k = 0; k < 2; k++)
			{
				//part->x[k] = part->x_half[k] + 0.5 * dt * part->v[k];
				//part->v[k] = part->v_half[k] + 0.5 * dt * part->a[k];
				part->x[k] += dt * part->v[k];
				part->v[k] += dt * part->a[k];
			}
		}
	}
	
	
	for (int k = 0; k < 2; k++) {
		part->a[k] = 0.0 + g[k];
	}
	//part->rho = part->rho_half + 0.5 * part->D * dt;
	part->rho += part->D * dt;
	part->P = B * (pow((part->rho / rho0), gamma) - 1);

	part->D = 0.0;
}


void SPH_main::reset_grid_count()
{
	for (auto& temp : grid_count) fill(temp.begin(), temp.end(), 0);

}

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


void SPH_main::store_initial(SPH_particle* part)
{
	part->rho_half = part->rho;
	part->v_half[0] = part->v[0];
	part->v_half[1] = part->v[1];
	part->x_half[0] = part->x[0];
	part->x_half[1] = part->x[1];

}

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


void SPH_main::time_dynamic()
{
	if (v_max != 0) dt_cfl = (h / v_max);
	if (a_max != 0) dt_a = sqrt(h / a_max);
	if (rho_max != 0) dt_f = h / (c0 * sqrt(pow((rho_max / rho0), gamma - 1)));

	if (dt_cfl <= dt_a && dt_cfl <= dt_f && dt_cfl > 1e-6)
	{
		dt = cfl * dt_cfl;
	}
	else if (dt_a <= dt_cfl && dt_a <= dt_f && dt_a > 1e-6)
	{
		dt = cfl * dt_a;
	}
	else if (dt_a <= dt_f && dt_a > 1e-6 && dt_cfl > 1e-6 )
	{
		dt = cfl * dt_a;
	}
	else if (dt_f > 1e-6)
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


double SPH_main::repulsion(SPH_particle *part, double &dist)
{
	double fd = pow((0.5 * dx / dist), 6) - 1.0;
	double temp_a = fd * 9.81 * 1.0 / dx;

    return temp_a;
}

