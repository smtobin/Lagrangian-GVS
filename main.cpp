#include "common.hpp"
#include "robot.hpp"

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
    
    ContinuumRobot robot(
        r,
        L,
        xi_star,
        E,
        nu,
        rho,
        strain_basis,
        n
    );

    robot.simulateRobot(200);
}