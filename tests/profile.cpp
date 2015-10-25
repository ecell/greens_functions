#include "../GreensFunction1DAbsAbs.hpp"
#include "../GreensFunction1DAbsSinkAbs.hpp"
#include "../GreensFunction1DRadAbs.hpp"
#include "../GreensFunction2DAbsSym.hpp"
#include "../GreensFunction2DRadAbs.hpp"
#include "../GreensFunction3D.hpp"
#include "../GreensFunction3DAbs.hpp"
#include "../GreensFunction3DSym.hpp"
#include "../GreensFunction3DAbsSym.hpp"
#include "../GreensFunction3DRadAbs.hpp"
#include "../GreensFunction3DRadInf.hpp"
#include <boost/random.hpp>
#include <iostream>
#include <time.h>

const int  loopnum(10000);

//using namespace greens_functions;

template <typename T_>
void test_drawTime(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1e0);
    Real bufT;

    clock_t start, end, dur;

    start = clock();
    for (int i=0; i < loopnum; i++)
    {
    bufT = gf.drawTime(randomarray[i]);
    }
    end = clock();

    std::cout << std::setw(30) << name << " ";
    std::cout << std::setw(15) << "drawtime: " << loopnum << " ";
    std::cout << "times: " << (double)(end - start) * 1000 / CLOCKS_PER_SEC << " msec" << std::endl;
//    std::cout << std::endl;
}

template <typename T_>
void test_drawR(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1e0);
    Real bufR;

    clock_t start, end, dur;

    start = clock();
    for (int i=0; i < loopnum; i++)
    {
      bufR = gf.drawR(randomarray[i], t);
    }
    end = clock();

    std::cout << std::setw(30) << name << " ";
    std::cout << std::setw(15) << "drawR: " << loopnum << " ";
    std::cout << "times: " << (double)(end - start) * 1000 / CLOCKS_PER_SEC << " msec" << std::endl;
//    std::cout << std::endl;
}
template <typename T_>
void test_theta(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1e0);
    Real const r(5e-1);
    Real bufTheta;
    clock_t start, end;

    start = clock();
    for (int i=0; i < loopnum; i++)
    {
    bufTheta = gf.drawTheta(randomarray[i], r, t);
    }
    end = clock();

    std::cout << std::setw(30) << name << " ";
    std::cout << std::setw(15) << "drawTheta: " << loopnum << " ";
    std::cout << "times: " << (double)(end - start) * 1000 / CLOCKS_PER_SEC << "msec" << std::endl;
//    std::cout << std::endl;
}

template <typename T_>
void test_event(T_ gf, std::string const name, Real* randomarray)
{
    Real const t(1e0);
    int eventType;
    clock_t start, end;

    start = clock();
    for (int i=0; i < loopnum; i++)
    {
    eventType = gf.drawEventType(randomarray[i], t);
    }
    end = clock();

    std::cout << std::setw(30) << name << " ";
    std::cout << std::setw(15) << "drawEvent: " << loopnum << " ";
    std::cout << "times: " << (double)(end - start) * 1000 / CLOCKS_PER_SEC << "msec" << std::endl;
//    std::cout << std::endl;
}

int main()
{
    Real const D(1e0);
    Real const a(1e1);
    Real const r0(1e0);
    Real const k(5e-1);
    Real const rsink(5e-1);
    Real const sigma(1e-1);
    Real randomarray[loopnum];


    boost::random::mt19937 rng(0);
    boost::random::uniform_real_distribution<Real> rand(0.0, 1.0);
    for(int i=0; i < loopnum; i++)
    {
    randomarray[i] = rand(rng);
    }

    test_drawTime(greens_functions::GreensFunction1DAbsAbs(D, r0, sigma, a), "GreensFunction1DAbsAbs", randomarray); //drawEventType
    test_drawTime(greens_functions::GreensFunction1DAbsSinkAbs(D, k, r0, rsink, sigma, a), "GreensFunction1DAbsSinkAbs", randomarray); //drawEventType
    test_drawTime(greens_functions::GreensFunction1DRadAbs(D, k, r0, sigma, a), "GreensFunction1DRadAbs", randomarray); //drawEventType
    test_drawTime(greens_functions::GreensFunction2DAbsSym(D, a), "GreensFunction2DAbsSym", randomarray); //
    test_drawTime(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray); //drawEventType, drawTheta
    test_drawTime(greens_functions::GreensFunction3D(D, a), "GreensFunction3D", randomarray); //drawTheta
    test_drawTime(greens_functions::GreensFunction3DAbs(D, r0, a), "GreensFunction3DAbs", randomarray); //drawEventType drawTheta
    test_drawTime(greens_functions::GreensFunction3DAbsSym(D, a), "GreensFunction3DAbsSym", randomarray); //
    test_drawTime(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray); //drawEventType, drawTheta
    test_drawTime(greens_functions::GreensFunction3DRadInf(D, k, r0, sigma), "GreensFunction3DRadInf", randomarray); //drawTheta
    std::cout << std::endl;

    test_drawR(greens_functions::GreensFunction1DAbsAbs(D, r0, sigma, a), "GreensFunction1DAbsAbs", randomarray); //drawEventType
    test_drawR(greens_functions::GreensFunction1DAbsSinkAbs(D, k, r0, rsink, sigma, a), "GreensFunction1DAbsSinkAbs", randomarray); //drawEventType
    test_drawR(greens_functions::GreensFunction1DRadAbs(D, k, r0, sigma, a), "GreensFunction1DRadAbs", randomarray); //drawEventType
    test_drawR(greens_functions::GreensFunction2DAbsSym(D, a), "GreensFunction2DAbsSym", randomarray); //
    test_drawR(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray); //drawEventType, drawTheta
    test_drawR(greens_functions::GreensFunction3D(D, a), "GreensFunction3D", randomarray); //drawTheta
    test_drawR(greens_functions::GreensFunction3DAbs(D, r0, a), "GreensFunction3DAbs", randomarray); //drawEventType drawTheta
    test_drawR(greens_functions::GreensFunction3DAbsSym(D, a), "GreensFunction3DAbsSym", randomarray); //
    test_drawR(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray); //drawEventType, drawTheta
    test_drawR(greens_functions::GreensFunction3DRadInf(D, k, r0, sigma), "GreensFunction3DRadInf", randomarray); //drawTheta
    std::cout << std::endl;

    test_theta(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray);
    test_theta(greens_functions::GreensFunction3D(D, a), "GreensFunction3D", randomarray);
    test_theta(greens_functions::GreensFunction3DAbs(D, r0, a), "GreensFunction3DAbs", randomarray);
    test_theta(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray);
    test_theta(greens_functions::GreensFunction3DRadInf(D, k, r0, sigma), "GreensFunction3DRadInf", randomarray);
    std::cout << std::endl;

    test_event(greens_functions::GreensFunction1DAbsAbs(D, r0, sigma, a), "GreensFunction1DAbsAbs", randomarray); //drawEventType
    test_event(greens_functions::GreensFunction1DAbsSinkAbs(D, k, r0, rsink, sigma, a), "GreensFunction1DAbsSinkAbs", randomarray); //drawEventType
    test_event(greens_functions::GreensFunction1DRadAbs(D, k, r0, sigma, a), "GreensFunction1DRadAbs", randomarray); //drawEventType
    test_event(greens_functions::GreensFunction2DRadAbs(D, k, r0, sigma, a), "GreensFunction2DRadAbs", randomarray); //drawEventType, drawTheta
//    test_event(greens_functions::GreensFunction3DAbs(D, r0, a), "GreensFunction3DAbs", randomarray); //drawEventType drawTheta
    test_event(greens_functions::GreensFunction3DRadAbs(D, k, r0, sigma, a), "GreensFunction3DRadAbs", randomarray); //drawEventType, drawTheta
    return 0;
}
