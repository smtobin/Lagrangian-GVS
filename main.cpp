#include "common.hpp"
#include "robot.hpp"
#include "visualization.hpp"

int main()
{
    // robot parameters
    Real r = 1e-3;
    Real L = 1;
    Vec6r xi_star(0, 0, 1, 1, 0, 0);
    Real E = 60e9;
    Real nu = 0.3;
    Real rho = 6500;

    // create strain basis
    unsigned n = 9;
    StrainBasisFunc strain_basis = [](Real s)
    {
        StrainBasisMat Phi(6, 9);
        Phi.col(0) = Vec6r(0, 0, 0, 1, 0, 0);
        Phi.col(1) = Vec6r(0, 0, 0, 0, 1, 0);
        Phi.col(2) = Vec6r(0, 0, 0, 0, 0, 1);
        Phi.col(3) = Vec6r(0, 0, 0, s, 0, 0);
        Phi.col(4) = Vec6r(0, 0, 0, 0, s, 0);
        Phi.col(5) = Vec6r(0, 0, 0, 0, 0, s);
        Phi.col(6) = Vec6r(0, 0, 0, s*s, 0, 0);
        Phi.col(7) = Vec6r(0, 0, 0, 0, s*s, 0);
        Phi.col(8) = Vec6r(0, 0, 0, 0, 0, s*s);
        return Phi;
    };
    
    ContinuumRobot robot1(
        r,
        L,
        xi_star,
        E,
        nu,
        rho,
        strain_basis,
        n
    );

    ContinuumRobot robot2 = robot1;

    Real h = 0.02;
    std::vector<VecXr> qs = robot1.simulateRobot(200, h, false, 1e-8);  // without C
    // robot2.simulateRobot(200, true, 1e-8);   // with C

    Visualization viz;
    viz.addGraphicsObject(&robot1, 30);
    viz.setup();
    viz.displayWindow();
    viz.interactorStart();
}