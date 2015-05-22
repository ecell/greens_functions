#include "GreensFunction1DAbsAbs.hpp"
#include "GreensFunction1DAbsSinkAbs.hpp"
#include "GreensFunction1DRadAbs.hpp"
#include "GreensFunction2DAbsSym.hpp"
#include "GreensFunction2DRadAbs.hpp"
#include "GreensFunction3D.hpp"
#include "GreensFunction3DAbs.hpp"
#include "GreensFunction3DSym.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "GreensFunction3DRadInf.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <time.h>

#define loopnum 1000

//using namespace greens_functions;

template <typename T_>
void test_func(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1.0);
    Real bufT, bufR;

    std::cout << name << std::endl;

    clock_t start, end, dur;
    start = clock();

    for (int i=0; i < 1000; i++)
    {
	bufT = gf.drawTime(randomarray[i]);
	bufR = gf.drawR(randomarray[i], t);
    }

    end = clock();

    std::cout << (double)(end - start) / 1000 << "msec" << std::endl;

//    std::cout << "drawTime: " << gf.drawTime(rand(rng)) << std::endl;
//    std::cout << "drawR: " << gf.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
}

template <typename T_>
void test_theta(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1.0);
    Real const r(5e-1);

    std::cout << name << std::endl;

    Real bufTheta;
    int eventType;

    clock_t start, end, dur;
    start = clock();

    for (int i=0; i < 1000; i++)
    {
	bufTheta = gf.drawTheta(randomarray[i], r, t);
	eventType = gf.drawEventType(randomarray[i], t);
    }

    end = clock();

    std::cout << (double)(end - start) / 1000 << "msec" << std::endl;

    std::cout << std::endl;
}

int main()
{
    Real const D(1.0);
    Real const a(1.e1);
    Real const r0(1.0);
    Real const k(0.5);
    Real const rsink(0.5);
    Real const sigma(1.e-1);
    
//    std::cout << CLOCKS_PER_SEC << std::endl;

    Real randomarray[loopnum];
    for(int i=0; i < loopnum; i++)
    {
	boost::random::mt19937 rng(0);
	boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);
	randomarray[i] = rand(rng);
    }

    test_func(greens_functions::GreensFunction1DAbsAbs(D, r0, sigma, a), "GreensFunction1DAbsAbs", randomarray); //drawEventType
    test_func(greens_functions::GreensFunction1DAbsSinkAbs(D, k, r0, rsink, sigma, a), "GreensFunction1DAbsSinkAbs", randomarray); //drawEventType
    test_func(greens_functions::GreensFunction1DRadAbs(D, k, r0, sigma, a), "GreensFunction1DRadAbs", randomarray); //drawEventType
    test_func(greens_functions::GreensFunction2DAbsSym(D, a), "GreensFunction2DAbsSym", randomarray); //
    test_func(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray); //drawEventType, drawTheta
    test_func(greens_functions::GreensFunction3D(D, a), "GreensFunction3D", randomarray); //drawTheta
    test_func(greens_functions::GreensFunction3DAbs(D, r0, a), "GreensFunction3DAbs", randomarray); //drawEventType drawTheta
    test_func(greens_functions::GreensFunction3DAbsSym(D, a), "GreensFunction3DAbsSym", randomarray); //
    test_func(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray); //drawEventType, drawTheta
    test_func(greens_functions::GreensFunction3DRadInf(D, k, r0, sigma), "GreensFunction3DRadInf", randomarray); //drawTheta

    test_theta(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray); //drawEventType, drawTheta
    test_theta(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray); //drawEventType, drawTheta

/*
    boost::random::mt19937 rng(0);
    boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);

    greens_functions::GreensFunction1DAbsAbs gf_1DAbsAbs(D, r0, sigma, a);
    greens_functions::GreensFunction1DAbsSinkAbs gf_1DAbsSinkAbs(D, k, r0, rsink, sigma, a);
    greens_functions::GreensFunction1DRadAbs gf_1DRadAbs(D, k, r0, sigma, a);
    greens_functions::GreensFunction2DAbsSym gf_2DAbsSym(D, a);
    greens_functions::GreensFunction2DRadAbs gf_2DRadAbs(D, k, r0, sigma, a);
    greens_functions::GreensFunction3D gf_3D(D, a);
    greens_functions::GreensFunction3DAbs gf_3DAbs(D, r0, a);
    greens_functions::GreensFunction3DSym gf_3DSym(D);
    greens_functions::GreensFunction3DAbsSym gf_3DAbsSym(D, a);
    greens_functions::GreensFunction3DRadAbs gf_3DRadAbs(D, k, r0, sigma, a);
    greens_functions::GreensFunction3DRadInf gf_3DRadInf(D, k, r0, sigma);

    Real const t(1.0);
    std::cout << "1DAbsAbs" << std::endl;
    std::cout << "drawTime: " << gf_1DAbsAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_1DAbsAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "1DAbsSinkAbs" << std::endl;
    std::cout << "drawTime: " << gf_1DAbsSinkAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_1DAbsSinkAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "1DRadAbs" << std::endl;
    std::cout << "drawTime: " << gf_1DRadAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_1DRadAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "2DAbsSym" << std::endl;
    std::cout << "drawTime: " << gf_2DAbsSym.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_2DAbsSym.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "2DRadAbs" << std::endl;
    std::cout << "drawTime: " << gf_2DRadAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_2DRadAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3D" << std::endl;
    std::cout << "drawTime: " << gf_3D.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3D.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3DAbs" << std::endl;
    std::cout << "drawTime: " << gf_3DAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3DAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3DSym" << std::endl;
    std::cout << "drawTime: " << gf_3DSym.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3DSym.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3DAbsSym" << std::endl;
    std::cout << "drawTime: " << gf_3DAbsSym.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3DAbsSym.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3DRadAbs" << std::endl;
    std::cout << "drawTime: " << gf_3DRadAbs.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3DRadAbs.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "3DRadInf" << std::endl;
    std::cout << "drawTime: " << gf_3DRadInf.drawTime(rand(rng)) << std::endl;
    std::cout << "drawR: " << gf_3DRadInf.drawR(rand(rng), t) << std::endl;
    std::cout << std::endl;
    std::cout << "end" << std::endl;
*/
    return 0;
}
