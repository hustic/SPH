#define BOOST_TEST_MODULE test SPH_2D
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/bind.hpp>

#include "SPH_2D.h"
#include "file_writer.h"

using namespace boost::unit_test;
namespace tt=boost::test_tools;

SPH_main domain;

int values[] = {};

BOOST_AUTO_TEST_SUITE(DomainConstructor);

BOOST_AUTO_TEST_CASE(TestDomainValues)
{
    domain.set_values();										//Set simulation parameters
    BOOST_TEST(true);
}

BOOST_AUTO_TEST_CASE(TestInitialiseGrid)
{
    domain.initialise_grid();									//initialise simulation grid
    BOOST_TEST(domain.search_grid.size() == domain.max_list[0]);
    for (int i = 0; i < domain.max_list[0]; i++)
        BOOST_TEST(domain.search_grid[i].size() == domain.max_list[1]);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(TestSplineFunctions);

BOOST_AUTO_TEST_CASE(TestCubicSplineLess)
{
    double dn[2] = { domain.h/2 , domain.h/2 };
    double q = sqrt(dn[0] * dn[0] + dn[1] * dn[1]) / domain.h;
    double tres = domain.cubic_spline(dn);
    double res = 10.0 * (1.0 - 1.5 * q * q + 0.75 * pow(q, 3)) / (7.0 * M_PI * domain.h * domain.h);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineGreater)
{
    double dn[2] = {domain.h, domain.h};
    double q = sqrt(dn[0] * dn[0] + dn[1] * dn[1]) / domain.h;
    double tres = domain.cubic_spline(dn);
    double res = 10  * 0.25 * pow((2 - q), 3) / (7 * M_PI * domain.h * domain.h);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineOutsite)
{
    double dn[2] = {2*domain.h, 2*domain.h};
    double tres = domain.cubic_spline(dn);
    double res = 0;
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstLess)
{
    double dn[2] = { domain.h/2 , domain.h/2 };
    double q = sqrt(dn[0] * dn[0] + dn[1] * dn[1]) / domain.h;
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = 10 * (-3 * q + 2.25 * q * q) / (7 * M_PI * pow(domain.h, 3));
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstGreater)
{
    double dn[2] = {domain.h, domain.h};
    double q = sqrt(dn[0] * dn[0] + dn[1] * dn[1]) / domain.h;
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = -10 * 0.75 * (2 - q) * (2 - q) / (7 * M_PI * pow(domain.h, 3));
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstOutside)
{
    double dn[2] = {2*domain.h, 2*domain.h};
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = 0;
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_SUITE_END();
