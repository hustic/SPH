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

/* 
 * Spline Functions Tests
 */
BOOST_AUTO_TEST_SUITE(TestSplineFunctions, * description("Testsing Cubic Spline Functions"));

BOOST_AUTO_TEST_CASE(TestCubicSplineLess,  * description("Distance between 0 and 1"))
{
    domain.h = 1.;
    double dn[2] = { 0.5 , 0. };
    double tres = domain.cubic_spline(dn);
    double res = 10.0 * (1.0 - 1.5 * 1/4 + 0.75 * 1/8) / (7.0 * M_PI);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineGreater, * description("Distance between 1 and 2"))
{
    domain.h = 1.;
    double dn[2] = {1, 0};
    double tres = domain.cubic_spline(dn);
    double res = 10  * 0.25 / (7 * M_PI);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineOutsite, * description("Distance outside the spline"))
{
    domain.h = 1.;
    double dn[2] = {2, 1};
    double tres = domain.cubic_spline(dn);
    double res = 0;
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstLess, \
* description("First derivative of distance between 0 and 1"))
{
    domain.h = 1.;
    double dn[2] = {0.5, 0};
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = 10 * (-3 * 0.5 + 2.25 * 1/4) / (7 * M_PI);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstGreater, \
* description("First derivative of distance between 1 and 2"))
{
    domain.h = 1.;
    double dn[2] = {1, 0};
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = -10 * 0.75 / (7 * M_PI);
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_CASE(TestCubicSplineFirstOutside, \
* description("First derivative of distance outside the spline"))
{
    domain.h = 1.;
    double dn[2] = {2., 1.};
    double tres = domain.cubic_spline_first_derivative(dn);
    double res = 0;
    BOOST_TEST(tres == res, tt::tolerance(0.0001));
}

BOOST_AUTO_TEST_SUITE_END();


/*
 * Gradient Function Tests
 */

BOOST_AUTO_TEST_SUITE(TestGradient, * description("Testing the Gradient Update Function"));

BOOST_AUTO_TEST_CASE(TestGradientLess, * description("Distance between 0 and 1"))
{
    SPH_particle* p1 = new SPH_particle();
    SPH_particle* p2 = new SPH_particle();
    
    domain.h = 1.;
    domain.mass = 1.;
    domain.mu = 1.;
    
    p1->v[0] = 1.;
    p1->v[1] = 0.;
    p1->rho = 1.;
    p1->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    p1->D = 0.;
    
    p2->v[0] = 0.;
    p2->v[1] = 0.;
    p2->rho = 1.;
    p2->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    
    
    double dn[2] = {0.5 , 0};
    
    domain.update_gradients(dn, p1, p2);
    
    double res = (10 * (-3 * 0.5 + 2.25 * 1/4) / (7 * M_PI));
    
    BOOST_TEST(2 * res == p1->a[0], tt::tolerance(0.0001));
    BOOST_TEST(res == p1->D, tt::tolerance(0.0001));
    BOOST_TEST(0 == p1->a[1], tt::tolerance(0.0001));
    
    std::cout << p1->D << " " << p1->a[0] << std::endl;
    
    std::cout << 2 * (10 * (-3 * 0.5 + 2.25 * 1/4) / (7 * M_PI)) << std::endl;
    
    delete p1;
    delete p2;
    
    // Need to check a and D
}

BOOST_AUTO_TEST_CASE(TestGradientGreater, * description("Distance between 1 and 2"))
{
    SPH_particle* p1 = new SPH_particle();
    SPH_particle* p2 = new SPH_particle();
    
    domain.h = 1.;
    domain.mass = 1.;
    domain.mu = 1.;
    
    p1->v[0] = 1.;
    p1->v[1] = 0.;
    p1->rho = 1.;
    p1->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    p1->D = 0.;
    
    p2->v[0] = 0.;
    p2->v[1] = 0.;
    p2->rho = 1.;
    p2->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    
    
    double dn[2] = {1 , 0};
    
    domain.update_gradients(dn, p1, p2);
    
    double res = -10 * 0.75 / (7 * M_PI);
    
    BOOST_TEST(2 * res == p1->a[0], tt::tolerance(0.0001));
    BOOST_TEST(res == p1->D, tt::tolerance(0.0001));
    BOOST_TEST(0 == p1->a[1], tt::tolerance(0.0001));
    
    std::cout << p1->D << " " << p1->a[0] << std::endl;
    
    std::cout << 2 * (-10 * 0.75 / (7 * M_PI)) << std::endl;
    
    
    
    // Need to check a and D
}

BOOST_AUTO_TEST_CASE(TestGradinetOutside, * description("Distance outside the spline"))
{
    SPH_particle* p1 = new SPH_particle();
    SPH_particle* p2 = new SPH_particle();
    
    domain.h = 1.;
    domain.mass = 1.;
    domain.mu = 1.;
    
    p1->v[0] = 1.;
    p1->v[1] = 0.;
    p1->rho = 1.;
    p1->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    p1->D = 0.;
    
    p2->v[0] = 0.;
    p2->v[1] = 0.;
    p2->rho = 1.;
    p2->P = 1.;
    p1->a[0] = 0.;
    p1->a[1] = 0.;
    
    
    double dn[2] = {0.5 , 0};
    
    domain.update_gradients(dn, p1, p2);
    
    double res = (10 * (-3 * 0.5 + 2.25 * 1/4) / (7 * M_PI));
    
    BOOST_TEST(2 * res == p1->a[0], tt::tolerance(0.0001));
    BOOST_TEST(res == p1->D, tt::tolerance(0.0001));
    BOOST_TEST(0 == p1->a[1], tt::tolerance(0.0001));
    
    std::cout << p1->D << " " << p1->a[0] << std::endl;
    
    std::cout << 2 * (10 * (-3 * 0.5 + 2.25 * 1/4) / (7 * M_PI)) << std::endl;
    
    
    
    // Need to check a and D
}

BOOST_AUTO_TEST_SUITE_END();
