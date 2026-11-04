#pragma once

#include "common.hpp"
#include "math.hpp"
#include "state.hpp"

#include <chrono>

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

struct ImplicitMidpointOutputChristoffel
{
    unsigned n;
    VecXr q;
    VecXr qdot;
    Real t;
    ContinuumRobotStateChristoffel tip_state;
    unsigned iter_to_convergence;
};

/////////////////////////////////////////////////////////////////////////////////////

struct ContinuumRobot
{
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
    ContinuumRobotStateChristoffel _state0_C;
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
        , _state0_C(n)
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
        {
            _s_ips[i] = Real(i) / Real(_s_ips.size() - 1) * _L;
        }

        // compute the stiffness matrix
        std::vector<MatXr> output = Math::RK4<MatXr>(MatXr::Zero(n,n), _s_ips, 
            [&](Real s, const MatXr& Kcur) {
                StrainBasisMat Phi = _strain_basis(s);
                MatXr Kp = Phi.transpose() * _Kd * Phi;
                return Kp;
            });
        _K = output.back();
    }

    std::vector<Mat4r> integratePoseAlongLength(const std::vector<Real>& s, const VecXr& q = _q)
    {
        return Math::RK4<Mat4r>(Mat4r::Identity(), s, [&](Real s, const Mat4r& Tcur) {
            StrainBasisMat Phi = _strain_basis(s);
            Vec6r xi = Phi*q + _xi_star;
            Mat4r Tp = Tcur * Math::se3_hat(xi);
            return Tp;
        });
    }

    ContinuumRobotState dynamicsODE(Real s, const ContinuumRobotState& state)
    {
        StrainBasisMat Phi = _strain_basis(s);
        Vec6r xi = Phi*_q + _xi_star;
        Vec6r xidot = Phi*_qdot;
        Mat6r adxi = Math::se3_ad(xi);
        Vec6r V = state.J*_qdot;
        Mat6r adV = Math::se3_ad(V);

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

    ContinuumRobotStateChristoffel dynamicsODEChristoffel(Real s, const ContinuumRobotStateChristoffel& state)
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
        StrainBasisMat Ap = adV*Phi - adxidot*state.J - adxi*state.A;
        MatXr Cp = 0.5*state.J.transpose() * (_Md*state.A - adV.transpose()*_Md*state.J + adbarMV*state.J);
        VecXr Qp = state.J.transpose() * wb;

        Vec3r p = state.T.block<3,1>(0,3);
        Vec6r xi_diff = xi - _xi_star;
        Real PEp = 0.5*xi_diff.transpose() * _Kd * xi_diff;
        PEp += _rho * _area * p.dot(_a_grav);
        Real KEp = 0.5*V.transpose()*_Md*V;

        ContinuumRobotStateChristoffel derivative {
            state.n,
            Tp,
            Jp,
            Ap,
            Mp,
            Cp,
            Qp,
            PEp,
            KEp
        };

        return derivative;
    }

    ImplicitMidpointOutput implicitMidpointStep(const VecXr& qk, const VecXr& qdk, Real tk, Real h, Real tol=1e-4)
    {
        /** TODO: (11/02/26) Make this function NOT change _q and _qdot. I think q and qd should be constant params given to dynamics ODE. */
        _q = qk + h/2*qdk;
        _qdot = qdk;

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

            Real qnorm = (q_new - _q).norm()/_n;
            Real qdnorm = (qd_new - _qdot).norm()/_n;

            _q = q_new;
            _qdot = qd_new;
            midpoint_output.iter_to_convergence = j;
            midpoint_output.tip_state = end;
            if (qnorm < tol && qdnorm < tol)
                break;
        }

        midpoint_output.n = _n;
        midpoint_output.q = 2*_q - qk;
        midpoint_output.qdot = 2*_qdot - qdk;
        midpoint_output.t = tk + h;

        return midpoint_output;
    }

    ImplicitMidpointOutputChristoffel implicitMidpointStepChristoffel(const VecXr& qk, const VecXr& qdk, Real tk, Real h, Real tol=1e-4)
    {
        /** TODO: (11/02/26) Make this function NOT change _q and _qdot. I think q and qd should be constant params given to dynamics ODE. */
        _q = qk + h/2*qdk;
        _qdot = qdk;

        ImplicitMidpointOutputChristoffel midpoint_output;

        unsigned MAX_ITER = 20;
        for (unsigned j = 0; j < MAX_ITER; j++)
        {
            // integrate ODE to get state along robot
            std::vector<ContinuumRobotStateChristoffel> output = 
                Math::RK4<ContinuumRobotStateChristoffel>(_state0_C, _s_ips, [this](Real s, const ContinuumRobotStateChristoffel& state) 
                {return this->dynamicsODEChristoffel(s, state);});

            const ContinuumRobotStateChristoffel& end = output.back();

            MatXr H = 2.0/h * end.M + 2*end.C + h/2.0*_K;
            VecXr b = 2.0/h * end.M * qdk + end.C*_qdot - _K*qk + end.Q;

            VecXr qd_new = H.partialPivLu().solve(b);
            VecXr q_new = qk + h/2.0*qd_new;

            Real qnorm = (q_new - _q).norm()/_n;
            Real qdnorm = (qd_new - _qdot).norm()/_n;

            _q = q_new;
            _qdot = qd_new;
            midpoint_output.iter_to_convergence = j;
            midpoint_output.tip_state = end;
            if (qnorm < tol && qdnorm < tol)
                break;
        }

        midpoint_output.n = _n;
        midpoint_output.q = 2*_q - qk;
        midpoint_output.qdot = 2*_qdot - qdk;
        midpoint_output.t = tk + h;

        return midpoint_output;
    }

    std::vector<VecXr> simulateRobot(unsigned num_steps, Real h, bool use_christoffel=false, Real tol=1e-4)
    {
        std::vector<VecXr> qs(num_steps+1);
        std::vector<Real> KEs(num_steps+1);
        std::vector<Real> PEs(num_steps+1);

        Real t = 0;

        VecXr qk = _q;
        VecXr qdk = _qdot;

        KEs[0] = 0;
        PEs[0] = 0;
        qs[0] = _q;

        auto t1 = std::chrono::high_resolution_clock::now();
        for (unsigned i = 0; i < num_steps; i++)
        {
            if (use_christoffel)
            {
                ImplicitMidpointOutputChristoffel output = implicitMidpointStepChristoffel(qk, qdk, t, h, tol);
                qk = output.q;
                qdk = output.qdot;
                t = output.t;

                KEs[i+1] = output.tip_state.KE;
                PEs[i+1] = output.tip_state.PE;
                qs[i+1] = qk;
            }
            else
            {
                ImplicitMidpointOutput output = implicitMidpointStep(qk, qdk, t, h, tol);
                qk = output.q;
                qdk = output.qdot;
                t = output.t;

                KEs[i+1] = output.tip_state.KE;
                PEs[i+1] = output.tip_state.PE;
                qs[i+1] = qk;
            }
            
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()/1.0e6;
        std::cout << "Simulating " << num_steps << " steps with tolerance " << tol << " took " << elapsed_ms << " ms." << std::endl;

        // for (unsigned i = 0; i < num_steps; i++)
        // {
        //     std::cout << "KE " << i << ": " << KEs[i] << std::endl;
        // }

        return qs;
    }



};