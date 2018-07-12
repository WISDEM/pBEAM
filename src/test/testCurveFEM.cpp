//
//  testBeam.cpp
//  pbeam
//
//  Created by Andrew Ning on 2/7/12.
//  Copyright (c) 2012 NREL. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <math.h>
#include "CurveFEM.h"

// Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
// Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
// pg. 168
TEST_CASE( "curvefem_fixed_beam_n_1" ){

    double E = 2.0;
    double I = 3.0;
    double L = 4.0;
    double A = 5.0;
    double rho = 6.0;
    double omegaRPM = 0.0;

    int n = 1;

    int nodes = n+1;

    Vector z(nodes);
    for (int i = 0; i < nodes; i++) {
        z(i) = L*i;
    }

    Vector EIx(nodes);
    for (int i = 0; i < nodes; i++) {
        EIx(i) = E * I;
    }

    Vector EIy(nodes);
    for (int i = 0; i < nodes; i++) {
        EIy(i) = E * I;
    }

    Vector EA(nodes);
    for (int i = 0; i < nodes; i++) {
        EA(i) = 1.0;
    }

    Vector GJ(nodes);
    for (int i = 0; i < nodes; i++) {
        GJ(i) = 1.0;
    }

    Vector rhoA(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoA(i) = rho * A;
    }

    Vector rhoJ(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoJ(i) = 1.0;
    }

    Vector theta(nodes);
    Vector precurv(nodes);
    Vector presweep(nodes);
    theta.setZero();
    precurv.setZero();
    presweep.setZero();
    
    CurveFEM mycurve = CurveFEM(omegaRPM, theta, z, precurv, presweep, rhoA, true);
    Vector freq = mycurve.frequencies(EA, EIx, EIy, GJ, rhoJ);

    Vector truth(nodes*6);
    truth << 0.012572866969753228, 0.012591740426479996, 0.015715395588976708, 0.015715395588976708, 0.02528529867098764, 0.02528529867098764, 0.06883002554331129, 0.06893334807998643, 0.09116741813672731, 0.09116741813672731, 0.1548391375450802, 0.1548391375450802;

    double tol_pct = 1e-7;
    for (int k=0; k<nodes*6; k++)
      REQUIRE(freq(k) == Approx(truth(k)).epsilon(tol_pct));

}


// Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
// Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
// pg. 168
TEST_CASE( "curvefem_fixed_beam_n_2" ){

    double E = 2.0;
    double I = 3.0;
    double L = 4.0;
    double A = 5.0;
    double rho = 6.0;
    double omegaRPM = 0.0;

    int n = 2;

    int nodes = n+1;

    Vector z(nodes);
    for (int i = 0; i < nodes; i++) {
        z(i) = L*i;
    }

    Vector EIx(nodes);
    for (int i = 0; i < nodes; i++) {
        EIx(i) = E * I;
    }

    Vector EIy(nodes);
    for (int i = 0; i < nodes; i++) {
        EIy(i) = E * I;
    }

    Vector EA(nodes);
    for (int i = 0; i < nodes; i++) {
        EA(i) = 1.0;
    }

    Vector GJ(nodes);
    for (int i = 0; i < nodes; i++) {
        GJ(i) = 1.0;
    }

    Vector rhoA(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoA(i) = rho * A;
    }

    Vector rhoJ(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoJ(i) = 1.0;
    }

    Vector theta(nodes);
    Vector precurv(nodes);
    Vector presweep(nodes);
    theta.setZero();
    precurv.setZero();
    presweep.setZero();
    
    CurveFEM mycurve = CurveFEM(omegaRPM, theta, z, precurv, presweep, rhoA, true);
    Vector freq = mycurve.frequencies(EA, EIx, EIy, GJ, rhoJ);

    Vector truth(nodes*6);
    truth << 0.003912151523474228, 0.003912151523474228, 0.005852976699216038, 0.012582302518525731, 0.02044674778626943, 0.02471287957376343, 0.02471287957376343, 0.025285550340850595, 0.025285550340850595, 0.032042058241816634, 0.06888168035399123, 0.0835842095344649, 0.0835842095344649, 0.09116764495898616, 0.09116764495898616, 0.11193550172703298, 0.24259762924185935, 0.24259762924185935;

    double tol_pct = 1e-7;
    for (int k=0; k<nodes*6; k++)
      REQUIRE(freq(k) == Approx(truth(k)).epsilon(tol_pct));
}



// Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
// Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
// pg. 168
TEST_CASE( "curvefem_fixed_beam_n_3" ){

    double E = 2.0;
    double I = 3.0;
    double L = 4.0;
    double A = 5.0;
    double rho = 6.0;
    double omegaRPM = 0.0;

    int n = 3;

    int nodes = n+1;

    Vector z(nodes);
    for (int i = 0; i < nodes; i++) {
        z(i) = L*i;
    }

    Vector EIx(nodes);
    for (int i = 0; i < nodes; i++) {
        EIx(i) = E * I;
    }

    Vector EIy(nodes);
    for (int i = 0; i < nodes; i++) {
        EIy(i) = E * I;
    }

    Vector EA(nodes);
    for (int i = 0; i < nodes; i++) {
        EA(i) = 1.0;
    }

    Vector GJ(nodes);
    for (int i = 0; i < nodes; i++) {
        GJ(i) = 1.0;
    }

    Vector rhoA(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoA(i) = rho * A;
    }

    Vector rhoJ(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoJ(i) = 1.0;
    }

    Vector theta(nodes);
    Vector precurv(nodes);
    Vector presweep(nodes);
    theta.setZero();
    precurv.setZero();
    presweep.setZero();
    
    CurveFEM mycurve = CurveFEM(omegaRPM, theta, z, precurv, presweep, rhoA, true);
    Vector freq = mycurve.frequencies(EA, EIx, EIy, GJ, rhoJ);

    Vector truth(nodes*6);
    truth << 0.0017380701632556908, 0.0017380701632562416, 0.003847212436706601, 0.010926964080897513, 0.010926964080897626, 0.012576854615234626, 0.012587751208210572, 0.02106152327278283, 0.02282612879817582, 0.02528527062288509, 0.02528527062288511, 0.03087567297329799, 0.030875672973298095, 0.06885185586578292, 0.06891150914730364, 0.06953079721389932, 0.06953079721389936, 0.09116748336205782, 0.09116748336205784, 0.12496139758839682, 0.1308572353442806, 0.13085723534428065, 0.2608788373591399, 0.2608788373591401;

    double tol_pct = 1e-7;
    for (int k=0; k<nodes*6; k++)
      REQUIRE(freq(k) == Approx(truth(k)).epsilon(tol_pct));

}

