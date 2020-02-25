#define BOOST_TEST_MODULE test SPH_2D
#include <boost/test/included/unit_test.hpp>

#include "SPH_2D.h"
#include "file_writer.h"

using namespace boost::unit_test;

SPH_main domain;

int values[] = {};


BOOST_AUTO_TEST_CASE(TestDomainValues)
{
    domain.set_values();										//Set simulation parameters
	
    BOOST_TEST(true);
}


BOOST_AUTO_TEST_CASE(Test_initialise_grid)
{
    domain.initialise_grid();									//initialise simulation grid
    BOOST_TEST(domain.search_grid.size() == domain.max_list[0]);
    for (int i = 0; i < domain.max_list[0]; i++)
        BOOST_TEST(domain.search_grid[i].size() == domain.max_list[1]);
}

BOOST_AUTO_TEST_CASE(TestStuff)
{
    domain.place_points(domain.min_x,domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

	domain.place_points(domain.min_x[0], domain.min_x[1], domain.max_x[0], domain.max_x[1]);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain


	domain.allocate_to_grid();									//needs to be called for each time step

	domain.neighbour_iterate(&domain.particle_list[100]);		//finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle
}
