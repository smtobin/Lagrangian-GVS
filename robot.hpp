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

    ContinuumRobotState() = default;

    ContinuumRobotState(unsigned n_)
    {
        n = n_;
        T = Mat4r::Identity();
        J = StrainBasisMat::Zero(6, n);
        a = Vec6r::Zero();
        M = MatXr::Zero(n, n);
        c = VecXr::Zero(n);
        Q = VecXr::Zero(n);
        PE = 0;
        KE = 0;
    }

    ContinuumRobotState(
        unsigned n_,
        const Mat4r& T_,
        const StrainBasisMat& J_,
        const Vec6r& a_,
        const MatXr& M_,
        const VecXr& c_,
        const VecXr& Q_,
        Real PE_,
        Real KE_)
        : n(n_), T(T_), J(J_), a(a_),
        M(M_), c(c_), Q(Q_), PE(PE_), KE(KE_)
    {}

    /** Overload operators for integration with RK4 */
    ContinuumRobotState operator+(const ContinuumRobotState& o) const
    {
        assert(n == o.n);
        return { n, T + o.T, J + o.J, a + o.a, M + o.M, c + o.c, Q + o.Q, PE + o.PE, KE + o.KE };
    }

    ContinuumRobotState operator-(const ContinuumRobotState& o) const
    {
        assert(n == o.n);
        return { n, T - o.T, J - o.J, a - o.a, M - o.M, c - o.c, Q - o.Q, PE - o.PE, KE - o.KE };
    }

    ContinuumRobotState operator*(Real s) const
    {
        return { n, s*T, s*J, s*a, s*M, s*c, s*Q, s*PE, s*KE };
    }

    ContinuumRobotState operator/(Real s) const
    {
        return { n, T/s, J/s, a/s, M/s, c/s, Q/s, PE/s, KE/s };
    }
};

ContinuumRobotState operator*(Real s, const ContinuumRobotState& state)
{
    return state*s;
}

/////////////////////////////////////////////////////////////////////////////////////

struct ImplicitMidpointOutput
{
    unsigned n;
    VecXr q;
    VecXr qdot;
    Real t;
    ContinuumRobotState tip_state;
    unsigned iter_to_convergence;
};

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

    /** Acceleration due to gravity */
    Vec3r _a_grav;

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

    /** Stiffness matrix - constant, integrated once at the beginning. */
    MatXr _K;

    /** Pre-computed things for integration along length of robot. */
    ContinuumRobotState _state0;
    std::vector<Real> _s_ips;

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
        , _a_grav(-9.81, 0, 0)
        , _strain_basis(strain_basis)
        , _n(n)
        , _state0(n)
        , _s_ips(30)
    {
        // initialize diagonal mass matrix
        _Md = Mat6r::Zero();
        _Md.diagonal() = Vec6r( _area*_rho, _area*_rho, _area*_rho, _rho*_Ixx, _rho*_Iyy, _rho*_Izz);

        _Kd = Mat6r::Zero();
        _Kd.diagonal() = Vec6r( _area*_G, _area*_G, _area*_E, _Ixx*_E, _Iyy*_E, _Izz*_G);

        // initialize generalized coordinates and their derivatives
        _q = VecXr::Zero(_n);
        _qdot = VecXr::Zero(_n);

        // initialize integration points for integration along robot length
        for (unsigned i = 0; i < _s_ips.size(); i++)
            _s_ips[i] = Real(i) / Real(_s_ips.size() - 1) * _L;

        // compute the stiffness matrix
        std::vector<MatXr> output = Math::RK4<MatXr>(MatXr::Zero(n,n), _s_ips, 
            [&](Real s, const MatXr& Kcur) {
                StrainBasisMat Phi = _strain_basis(s);
                return Phi.transpose() * _Kd * Phi;
            });
        _K = output.back();
        
    }

    ContinuumRobotState dynamicsODE(Real s, const ContinuumRobotState& state)
    {
        StrainBasisMat Phi = _strain_basis(s);
        Vec6r xi = Phi*_q + _xi_star;
        Vec6r xidot = Phi*_qdot;
        Mat6r adxi = Math::se3_ad(xi);
        Mat6r adxidot = Math::se3_ad(xidot);
        Vec6r V = state.J*_qdot;
        Mat6r adV = Math::se3_ad(V);
        Mat6r adbarMV = Math::se3_adbar(_Md*V);

        Mat3r R = state.T.block<3,3>(0,0);
        Vec6r wb = Vec6r::Zero();
        wb.head<3>() = R.transpose() * _a_grav * _rho * _area;

        Mat4r Tp = state.T * Math::se3_hat(xi);
        StrainBasisMat Jp = Phi - adxi*state.J;
        MatXr Mp = state.J.transpose() * _Md * state.J;
        Vec6r ap = adV*xidot - adxi*state.a;
        VecXr cp = state.J.transpose() * (_Md*state.a - adV.transpose() * _Md * V);
        VecXr Qp = state.J.transpose() * wb;

        Vec3r p = state.T.block<3,1>(0,3);
        Vec6r xi_diff = xi - _xi_star;
        Real PEp = 0.5*xi_diff.transpose() * _Kd * xi_diff;
        PEp += _rho * _area * p.dot(_a_grav);
        Real KEp = 0.5*V.transpose()*_Md*V;

        ContinuumRobotState derivative {
            state.n,
            Tp,
            Jp,
            ap,
            Mp,
            cp,
            Qp,
            PEp,
            KEp
        };

        return derivative;
    }

    ImplicitMidpointOutput implicitMidpointStep(const VecXr& qk, const VecXr& qdk, Real tk, Real h, Real tol=1e-4)
    {
        VecXr q = qk + h/2*qdk;
        VecXr qd = qdk;

        ImplicitMidpointOutput midpoint_output;

        unsigned MAX_ITER = 20;
        for (unsigned j = 0; j < MAX_ITER; j++)
        {
            // integrate ODE to get state along robot
            std::vector<ContinuumRobotState> output = 
                Math::RK4<ContinuumRobotState>(_state0, _s_ips, [this](Real s, const ContinuumRobotState& state) {return this->dynamicsODE(s, state);});
            const ContinuumRobotState& end = output.back();

            MatXr H = 2.0/h * end.M + h/2.0*_K;
            VecXr b = 2.0/h * end.M * qdk - end.c - _K*qk + end.Q;

            VecXr qd_new = H.partialPivLu().solve(b);
            VecXr q_new = qk + h/2.0*qd_new;

            Real qnorm = (q_new - q).norm()/_n;
            Real qdnorm = (qd_new - qd).norm()/_n;

            q = q_new;
            qd = qd_new;
            midpoint_output.iter_to_convergence = j;
            midpoint_output.tip_state = end;
            if (qnorm < 1e-4 && qdnorm < 1e-4)
                break;
        }

        midpoint_output.n = _n;
        midpoint_output.q = 2*q - qk;
        midpoint_output.qdot = 2*qd - qdk;
        midpoint_output.t = tk + h;

        return midpoint_output;
    }

    void simulateRobot(unsigned num_steps)
    {
        std::vector<Real> KEs(num_steps);
        std::vector<Real> PEs(num_steps);

        Real t = 0;
        Real h = 0.02;
        for (unsigned i = 0; i < num_steps; i++)
        {
            ImplicitMidpointOutput output = implicitMidpointStep(_q, _qdot, t, h);
            _q = output.q;
            _qdot = output.qdot;
            t = output.t;

            KEs[i] = output.tip_state.KE;
            PEs[i] = output.tip_state.PE;
        }

        for (unsigned i = 0; i < num_steps; i++)
        {
            std::cout << "KE " << i << ": " << KEs[i] << std::endl;
        }
    }



};