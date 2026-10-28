#pragma once

#include "common.hpp"
#include "math.hpp"

#include <functional>

struct ContinuumRobotState
{
    /** Number of generalized coordinates */
    unsigned n;
    /** Pose */
    Mat4r T;
    /** Pose Jacobian (6xn) */
    StrainBasisMat J;
    /** Jdot * qdot (?) */
    Vec6r a;
    /** Mass matrix (nxn) */
    MatXr M;
    /** Coriolis terms (C*qdot), (nx1) */
    VecXr c;
    /** Generalized force (nx1) */
    VecXr Q;
    /** Potential energy */
    Real PE;
    /** Kinetic energy */
    Real KE;

    /** Overload operators for integration with RK4 */
    ContinuumRobotState operator+(const ContinuumRobotState& o)
    {
        assert(n == o.n);
        return { n, T + o.T, J + o.J, a + o.a, M + o.M, c + o.c, Q + o.Q, PE + o.PE, KE + o.KE };
    }

    ContinuumRobotState operator-(const ContinuumRobotState& o)
    {
        assert(n == o.n);
        return { n, T - o.T, J - o.J, a - o.a, M - o.M, c - o.c, Q - o.Q, PE - o.PE, KE - o.KE };
    }

    ContinuumRobotState operator*(Real s)
    {
        return { n, s*T, s*J, s*a, s*M, s*c, s*Q, s*PE, s*KE };
    }

    ContinuumRobotState operator/(Real s)
    {
        return { n, T/s, J/s, a/s, M/s, c/s, Q/s, PE/s, KE/s };
    }
};

ContinuumRobotState operator*(Real s, const ContinuumRobotState& state)
{
    return state*s;
}

/////////////////////////////////////////////////////////////////////////////////////

class ContinuumRobot
{
private:
    /** Geometric properties */
    Real _r;        // radius
    Real _L;        // length
    Real _area;     // cross-sectional area
    Real _Ixx;      // second moment of area, x component
    Real _Iyy;      // second moment of area, y component
    Real _Izz;      // second moment of area, z component
    Vec6r _xi_star; // rest-state strain

    /** Material properties */
    Real _E;
    Real _nu;
    Real _G;
    Real _rho;

    /** Digaonal (body frame) mass matrix */
    Mat6r _Md;
    
    /** Diagonal (body frame) stiffness matrix */
    Mat6r _Kd;

    /** Represents the strain basis functions
     * Can be evaluated at a given s to obtain the Phi matrix
     */
    StrainBasisFunc _strain_basis;

    /** Number of states */
    unsigned _n;

    /** Generalized coordinates and their time derivatives */
    VecXr _q;
    VecXr _qdot;

public:
    ContinuumRobot(Real r, Real L, const Vec6r& xi_star, Real E, Real nu, Real rho, StrainBasisFunc strain_basis, unsigned n)
        : _r(r)
        , _L(L)
        , _area(M_PI*r*r)
        , _Ixx(M_PI*r*r*r*r/4.0)
        , _Iyy(_Ixx)
        , _Izz(_Ixx)
        , _xi_star(xi_star)
        , _E(E)
        , _nu(nu)
        , _G(E/(2*(1+nu)))
        , _rho(rho)
        , _strain_basis(strain_basis)
        , _n(n)
    {
        // initialize diagonal mass matrix
        _Md = Mat6r::Zero();
        _Md.diagonal() = Vec6r( _area*_rho, _area*_rho, _area*_rho, _rho*_Ixx, _rho*_Iyy, _rho*_Izz);

        _Kd = Mat6r::Zero();
        _Kd.diagonal() = Vec6r( _area*_G, _area*_G, _area*_E, _Ixx*_E, _Iyy*_E, _Izz*_G);

        // initialize generalized coordinates and their derivatives
        _q = VecXr::Zero(_n);
        _qdot = VecXr::Zero(_n);
    }

    ContinuumRobotState dynamicsODE(Real s, const ContinuumRobotState& state)
    {
        StrainBasisMat Phi = _strain_basis(s);
        Vec6r xi = Phi*_q + _xi_star;
        Vec6r xidot = Phi*_qdot;
        Mat6r adxi = Math::se3_ad(xi);
        Mat6r adxidot = Math::se3_ad(xidot);
        Vec6r V = J*_qdot;
        Mat6r adV = Math::se3_ad(V);
        Mat6r adbarMV = Math::se3_adbar(_Md*V);

        /** TODO: finish this function */
    }


};