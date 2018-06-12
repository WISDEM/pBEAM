import numpy as np
import numpy.testing as npt
import unittest
import _pBEAM as pb

class TestPyBeam(unittest.TestCase):
    
    def testCantileverDeflection(self):
        # Test data from "Finite Element Structural Analysis", Yang, pg. 145
        E = 2.0
        I = 3.0
        L = 4.0
        p0 = 5.0

        nodes = 2

        Px = np.array([-p0, 0.0])
        Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData()

        base = pb.BaseData(np.ones(6), 1.0)

        z = np.array([0.0, L])
        EIx = EIy = E*I*np.ones(nodes)
        EA = E*np.ones(nodes)
        GJ = rhoA = rhoJ = np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        (dx, dy, dz, dtheta_x, dtheta_y, dtheta_z) = beam.displacement()

        self.assertAlmostEqual(dx[0], 0.0, 8)
        self.assertAlmostEqual(dx[1], -p0 * L**3.0 / E / I * L / 30.0, 8)
        self.assertAlmostEqual(dtheta_y[0], 0.0, 8)
        self.assertAlmostEqual(dtheta_y[1],-p0 * L**3.0 / E / I * 1 / 24.0, 8)

    def testTaperedDeflections(self):
        # Test data from "Finite Element Structural Analysis", Yang, pg. 180
        E = 2.0
        I = 3.0
        L = 4.0
        P = 5.0

        nodes = 3

        Px = Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData(0.0, np.zeros(3), np.zeros(6), [-P, 0.0, 0.0], np.zeros(3))

        base = pb.BaseData(np.ones(6), 1.0)

        z = np.array([0.0, 0.5*L, L])
        EIx = EIy = E*I*np.array([9.0, 5.0, 1.0])
        EA = E*np.ones(nodes)
        GJ = rhoA = rhoJ = np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        (dx, dy, dz, dtheta_x, dtheta_y, dtheta_z) = beam.displacement()

        tol = 1e-8
        tol_pct_1 = 0.17;
        tol_pct_2 = 0.77;
        self.assertAlmostEqual(dx[0], 0.0, 8)
        self.assertAlmostEqual(dx[-1], -0.051166*P*L**3/E/I, delta=tol_pct_1)
        self.assertAlmostEqual(dtheta_y[0], 0.0, 8)
        self.assertAlmostEqual(dtheta_y[-1], -0.090668*P*L**2/E/I, delta=tol_pct_2)

    def testFreqFree_FreeBeam_n1(self):
        # Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
        # Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
        # pg. 168
        E = 2.0
        I = 3.0
        L = 4.0
        A = 5.0
        rho = 6.0

        n = 1
        nodes = n+1
        
        Px = Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData()

        base = pb.BaseData()
        
        z = L * np.arange(nodes)
        EIx = EIy = E*I*np.ones(nodes)
        EA = GJ = rhoJ = np.ones(nodes)
        rhoA = rho*A*np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        nFreq = 100
        freq = beam.naturalFrequencies(nFreq)
        
        m = rho * A
        alpha = m * (n*L)**4.0 / (840.0 * E * I)

        tol_pct = 5e-6 * 100
        expect = np.sqrt(0.85714 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[1], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[2], expect, delta=tol_pct)
        expect = np.sqrt(10.0 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[4], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[5], expect, delta=tol_pct)



    def testFreqFree_FreeBeam_n2(self):
        # Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
        # Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
        # pg. 168
        E = 2.0
        I = 3.0
        L = 4.0
        A = 5.0
        rho = 6.0

        n = 2
        nodes = n+1
        
        Px = Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData()

        base = pb.BaseData()
        
        z = L * np.arange(nodes)
        EIx = EIy = E*I*np.ones(nodes)
        EA = GJ = rhoJ = np.ones(nodes)
        rhoA = rho*A*np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        nFreq = 100
        freq = beam.naturalFrequencies(nFreq)
        
        m = rho * A
        alpha = m * (n*L)**4.0 / (840.0 * E * I)

        tol_pct = 6e-6 * 100
        expect = np.sqrt(0.59858 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[1], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[2], expect, delta=tol_pct)
        expect = np.sqrt(5.8629 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[5], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[6], expect, delta=tol_pct)
        expect = np.sqrt(36.659 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[8], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[9], expect, delta=tol_pct)
        expect = np.sqrt(93.566 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[10], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[11], expect, delta=tol_pct)
        

    def testFreqFree_FreeBeam_n3(self):
        # Test data from "Consistent Mass Matrix for Distributed Mass Systmes", John Archer,
        # Journal of the Structural Division Proceedings of the American Society of Civil Engineers,
        # pg. 168
        E = 2.0
        I = 3.0
        L = 4.0
        A = 5.0
        rho = 6.0

        n = 3
        nodes = n+1
        
        Px = Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData()

        base = pb.BaseData()
        
        z = L * np.arange(nodes)
        EIx = EIy = E*I*np.ones(nodes)
        EA = GJ = rhoJ = np.ones(nodes)
        rhoA = rho*A*np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        nFreq = 100
        freq = beam.naturalFrequencies(nFreq)
        
        m = rho * A
        alpha = m * (n*L)**4.0 / (840.0 * E * I)

        tol_pct = 6e-6 * 100
        expect = np.sqrt(0.59919 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[1], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[2], expect, delta=tol_pct)
        expect = np.sqrt(4.5750 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[5], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[6], expect, delta=tol_pct)
        expect = np.sqrt(22.010 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[8], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[9], expect, delta=tol_pct)
        expect = np.sqrt(70.920 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[11], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[12], expect, delta=tol_pct)
        expect = np.sqrt(265.91 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[14], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[15], expect, delta=tol_pct)
        expect = np.sqrt(402.40 / alpha) / (2*np.pi)
        self.assertAlmostEqual(freq[16], expect, delta=tol_pct)
        self.assertAlmostEqual(freq[17], expect, delta=tol_pct)

    def testBucklingEuler(self):
        # unit test data from Euler's buckling formula for a clamped/free beam

        E = 2.0
        I = 3.0
        L = 4.0

        n = 3
        nodes = n+1
        
        Px = Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        tip = pb.TipData()

        base = pb.BaseData(np.ones(6), 1.0)
        
        z = L * np.arange(nodes)
        EIx = EIy = E*I*np.ones(nodes)
        EA = GJ = rhoJ = rhoA = np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        beam = pb.Beam(sec, loads, tip, base)

        (Pcr_x, Pcr_y) = beam.criticalBucklingLoads()

        Ltotal = n*L
        tol_pct = 0.011
        expect = E * I * (0.5*np.pi/Ltotal)**2.0
        self.assertAlmostEqual(Pcr_x, expect, 2)
        self.assertAlmostEqual(Pcr_y, expect, 2)

    def testShearBendingSimple(self):
        # Test data from "Mechanical of Materials", Gere, 6th ed., pg. 273
        # cantilevered beam with linear distributed load

        L = 10.0
        q0 = 3.0

        n = 1
        nodes = n+1

        tip = pb.TipData()

        base = pb.BaseData(np.ones(6), 1.0)
        
        z = np.arange(nodes, dtype=np.float64) * (L / n)
        EIx = EIy = EA = GJ = rhoJ = rhoA = np.ones(nodes)
        sec = pb.SectionData(nodes, z, EA, EIx, EIy, GJ, rhoA, rhoJ)

        Px = q0*(1 - z/L)
        Py = Pz = np.zeros(nodes)
        loads = pb.Loads(nodes, Px, Py, Pz)

        beam = pb.Beam(sec, loads, tip, base)

        (Vx, Vy, Fz, Mx, My, Tz) = beam.shearAndBending()

        tol_pct = 1e-8
        Vx_expect = np.polyval([q0*L/2.0, -q0*L, q0*L/2.0], [0.0, 1.0])
        My_expect = np.polyval([-q0*L*L/6.0, 3.0*q0*L*L/6.0, -3.0*q0*L*L/6.0, q0*L*L/6.0], [0.0, 1.0])
        npt.assert_almost_equal(Vx, Vx_expect)
        npt.assert_almost_equal(My, My_expect)

        
        
'''    

// Test data from "Mechanical of Materials", Gere, 6th ed., pg. 288
// cantilevered beam with two point loads
TEST_CASE( "shear_bending_simple_pt" ){

    double L = 10.0;
    double P1 = 2.0;
    double P2 = 3.0;

    int n = 3;

    int nodes = n+1;

    Vector z(nodes);
    for (int i = 0; i < nodes; i++) {
        z(i) = (double)i/(nodes-1) * L;
    }

    Vector EIx(nodes);
    for (int i = 0; i < nodes; i++) {
        EIx(i) = 1.0;
    }

    Vector EIy(nodes);
    for (int i = 0; i < nodes; i++) {
        EIy(i) = 1.0;
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
        rhoA(i) = 1.0;
    }

    Vector rhoJ(nodes);
    for (int i = 0; i < nodes; i++) {
        rhoJ(i) = 1.0;
    }

    Vector ESx(nodes), ESy(nodes), EIxy(nodes);
    ESx.setZero();
    ESy.setZero();
    EIxy.setZero();

    Vector Px(nodes);
    Vector Py(nodes);
    Vector Pz(nodes);
    Px.setZero();
    Py.setZero();
    Pz.setZero();

    Vector Fx_pt(nodes), Fy_pt(nodes), Fz_pt(nodes), Mx_pt(nodes), My_pt(nodes), Mz_pt(nodes);
    Fx_pt.setZero();
    Fy_pt.setZero();
    Fz_pt.setZero();
    Mx_pt.setZero();
    My_pt.setZero();
    Mz_pt.setZero();

    Fx_pt(1) = -P2;
    Fx_pt(3) = -P1;

    Loads loads = Loads(Px, Py, Pz, Fx_pt, Fy_pt, Fz_pt, Mx_pt, My_pt, Mz_pt);

    TipData tip;

    BaseData base;
    for (int i = 0; i < 6; i++) {
        base.rigid[i] = true;
    }

    SectionData sec = SectionData(z, EA, EIx, EIy, GJ, rhoA, rhoJ);

    Beam beam = Beam(sec, loads, tip, base);

    PolyVec Vx, Vy, Fz, Mx, My, Tz;

    beam.shearAndBending(Vx, Vy, Fz, Mx, My, Tz);


    double tol_pct = 1e-8;
    self.assertAlmostEqual(Vx[0].length()== 3);
    self.assertAlmostEqual(Vx[0](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[0](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[0](2) == Approx(-P1-P2).epsilon(tol_pct));

    self.assertAlmostEqual(Vx[1](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[1](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[1](2) == Approx(-P1).epsilon(tol_pct));

    self.assertAlmostEqual(Vx[2](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[2](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Vx[2](2) == Approx(-P1).epsilon(tol_pct));

    double b = L/3.0;
    double a = 2.0/3.0*L;

    self.assertAlmostEqual(Mx[0](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[0](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[0](2) == Approx(-P1*a + P1*L + P2*b).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[0](3) == Approx(-P1*L - P2*b).epsilon(tol_pct));

    self.assertAlmostEqual(Mx[1](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[1](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[1](2) == Approx(-0.5*P1*a + P1*a).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[1](3) == Approx(-P1*a).epsilon(tol_pct));

    self.assertAlmostEqual(Mx[2](0) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[2](1) == Approx(0.0).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[2](2) == Approx(0.5*P1*a).epsilon(tol_pct));
    self.assertAlmostEqual(Mx[2](3) == Approx(-0.5*P1*a).epsilon(tol_pct));


}


// Test data from Rick Damiani's ANSYS model
// tower with springs at base, and offset mass
TEST_CASE( "ricks_tower" ){

    double dmp = 6.2;
    double dtw_base = 6.424;
    double dtw_top = dtw_base * 0.55;

    double E = 210.0e9;
    double G = 80.769e9;
    double rho = 8502.0;

    double E_grout = 3.9e10;
    double rho_grout = 2500.0;
    double G_grout = 1.5e10;

    Vector z(100);
    z(0) = 0.0;
    PolyVec EIx(100);
    PolyVec EIy(100);
    PolyVec EA(100);
    PolyVec GJ(100);
    PolyVec rhoA(100);
    PolyVec rhoJ(100);



    int idx = 0;

    // monopile bottom
    int n_mp_bot = 2;
    double I = np.pi / 8.0 * pow(dmp, 3) * dmp/76.938;
    double A = np.pi * dmp * dmp/76.938;
    for (int i = 1; i <= n_mp_bot; i++){
        z(idx+i) = z(idx) + (double)i/n_mp_bot*(2.0*dmp);
        EIx[idx+i-1] = Poly(1, E*I);
        EIy[idx+i-1] = Poly(1, E*I);
        EA[idx+i-1] = Poly(1, E*A);
        GJ[idx+i-1] = Poly(1, G*(I+I));
        rhoA[idx+i-1] = Poly(1, rho*A);
        rhoJ[idx+i-1] = Poly(1, rho*(I+I));
    }
    idx += n_mp_bot;

    // monopile bottom
    int n_mp_top = 2;
    I = np.pi / 8.0 * pow(dmp, 3) * dmp/100.0;
    A = np.pi * dmp * dmp/100.0;
    for (int i = 1; i <= n_mp_top; i++){
        z(idx+i) = z(idx) + (double)i/n_mp_top*(25.0 + 5.0 - 0.5*dmp - 2.0*dmp);
        EIx[idx+i-1] = Poly(1, E*I);
        EIy[idx+i-1] = Poly(1, E*I);
        EA[idx+i-1] = Poly(1, E*A);
        GJ[idx+i-1] = Poly(1, G*(I+I));
        rhoA[idx+i-1] = Poly(1, rho*A);
        rhoJ[idx+i-1] = Poly(1, rho*(I+I));
    }
    idx += n_mp_top;

    // overlap
    int n_ov = 1;
    double I1 = np.pi / 8.0 * pow(dtw_base, 3) * 0.062;
    double I2 = np.pi / 8.0 * pow(dmp, 3) * dmp/100.0;
    double I3 = np.pi / 64.0 * (pow(dtw_base-0.062, 4) - pow(dmp+dmp/100.0, 4));
    double A1 = np.pi * dtw_base * 0.062;
    double A2 = np.pi * dmp * dmp/100.0;
    double A3 = np.pi / 4.0 * (pow(dtw_base-0.062, 2) - pow(dmp+dmp/100.0, 2));
    for (int i = 1; i <= n_ov; i++){
        z(idx+i) = z(idx) + (double)i/n_ov*(0.5*dmp);
        EIx[idx+i-1] = Poly(1, E*I1 + E*I2 + E_grout*I3);
        EIy[idx+i-1] = Poly(1, E*I1 + E*I2 + E_grout*I3);
        EA[idx+i-1] = Poly(1, E*A1 + E*A2 + E_grout*A3);
        GJ[idx+i-1] = Poly(1, G*2*I1 + G*2*I2 + G_grout*2*I3);
        rhoA[idx+i-1] = Poly(1, rho*A1 + rho*A2 + rho_grout*A3);
        rhoJ[idx+i-1] = Poly(1, rho*2*I1 + rho*2*I2 + rho_grout*2*I3);
    }
    idx += n_ov;

    // transition
    int n_ts = 4;
    I = np.pi / 8.0 * pow(dtw_base, 3) * 0.062;
    A = np.pi * dtw_base * 0.062;
    for (int i = 1; i <= n_ts; i++){
        z(idx+i) = z(idx) + (double)i/n_ts*(18.05 - 0.5*dmp);
        EIx[idx+i-1] = Poly(1, E*I);
        EIy[idx+i-1] = Poly(1, E*I);
        EA[idx+i-1] = Poly(1, E*A);
        GJ[idx+i-1] = Poly(1, G*(I+I));
        rhoA[idx+i-1] = Poly(1, rho*A);
        rhoJ[idx+i-1] = Poly(1, rho*(I+I));
    }
    idx += n_ts;

    // tower
    int n_tw = 16;
    for (int i = 1; i <= n_tw; i++){
        z(idx+i) = z(idx) + (double)i/n_tw*(88.37);
        double d_bot = dtw_base + (double)(i-1)/n_tw*(dtw_top-dtw_base);
        double d_top = dtw_base + (double)(i)/n_tw*(dtw_top-dtw_base);
        Poly d = Poly(2, d_top-d_bot, d_bot);
        Poly t = d / 120.0;
        Poly Itw = np.pi / 8.0 * d * d * d* t;
        Poly Atw = np.pi * d * t;
        EIx[idx+i-1] = E * Itw;
        EIy[idx+i-1] = E * Itw;
        EA[idx+i-1] = E * Atw;
        GJ[idx+i-1] = G * 2*Itw;
        rhoA[idx+i-1] = rho * Atw;
        rhoJ[idx+i-1] = rho * 2*Itw;
    }
    idx += n_tw;

    int nodes = n_mp_bot + n_mp_top + n_ov + n_ts + n_tw + 1;
    z.resize(nodes);
    EIx.resize(nodes-1);
    EIy.resize(nodes-1);
    EA.resize(nodes-1);
    GJ.resize(nodes-1);
    rhoA.resize(nodes-1);
    rhoJ.resize(nodes-1);

    PolyVec Px(nodes-1);
    PolyVec Py(nodes-1);
    PolyVec Pz(nodes-1);
    //Px.setZero();
    //Py.setZero();
    //Pz.setZero();

    Vector Fx_pt(nodes), Fy_pt(nodes), Fz_pt(nodes), Mx_pt(nodes), My_pt(nodes), Mz_pt(nodes);
    Fx_pt.setZero();
    Fy_pt.setZero();
    Fz_pt.setZero();
    Mx_pt.setZero();
    My_pt.setZero();
    Mz_pt.setZero();

    double m = 5.7380e5;
    double Ixx = 86.579e6;
    double Iyy = 53.530e6;
    double Izz = 58.112e6;
    double Itip[] = {Ixx, Iyy, Izz, 0.0, 0.0, 0.0};
    double cm[] = {0.0, 0.0, 2.34};
    double F[] = {2000.0e3, 0.0, 0.0};
    double M[] = {0.0, 0.0, 0.0};

//    m = 0.0;
//    Itip[0] = 0.0;
//    Itip[1] = 0.0;
//    Itip[2] = 0.0;

    TipData tip(m, cm, Itip, F, M);

    double kx = 4.72e8;
    double ktx = 1.27e11;
    BaseData base(kx, ktx, kx, ktx, 999, 999, 999);
    //BaseData base(999, 999, 999, 999, 999, 999, 999);

    PolynomialSectionData sec = PolynomialSectionData(z, EA, EIx, EIy, GJ, rhoA, rhoJ);
    PolynomialLoads loads = PolynomialLoads(Px, Py, Pz, Fx_pt, Fy_pt, Fz_pt, Mx_pt, My_pt, Mz_pt);
    Beam beam = Beam(sec, loads, tip, base);

    int nFreq = 100;
    Vector freq(nFreq);

//    beam.computeNaturalFrequencies(nFreq, freq);
//    double mass = beam.computeMass();
//    double Pcr_x, Pcr_y;
//    beam.computeMinBucklingLoads(Pcr_x, Pcr_y);

//    std::cout << mass << std::endl;
//    for (int i = 0; i < 20; i++){
//        std::cout << freq(i) << std::endl;
//    }
//    std::cout << Pcr_x << std::endl;



//    beam.computeAxialStress(<#Vector &x#>, <#Vector &y#>, <#Vector &z#>, <#Vector &E#>, <#Vector &sigma_axial#>)

//    double tol_pct = 6e-6*100;
    //self.assertAlmostEqual(freq(1), sqrt(0.59919 / alpha) / (2*np.pi)).epsilon(tol_pct));

}
'''        
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestPyBeam))
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner().run(suite())
