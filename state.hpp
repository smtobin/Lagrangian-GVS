#pragma once

#include "common.hpp"
#include "math.hpp"

#include <functional>
#include <iostream>

/** Integration state along length of robot, without Christoffel-consistent C matrix. */
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

std::ostream& operator<<(std::ostream& os, const ContinuumRobotState& state)
{
    os << "ContinuumRobotState {\n"
       << "  n  = " << state.n << "\n"
       << "  T  =\n" << state.T << "\n"
       << "  J  =\n" << state.J << "\n"
       << "  a  =\n" << state.a.transpose() << "\n"
       << "  M  =\n" << state.M << "\n"
       << "  c  =\n" << state.c.transpose() << "\n"
       << "  Q  =\n" << state.Q.transpose() << "\n"
       << "  PE = " << state.PE << "\n"
       << "  KE = " << state.KE << "\n"
       << "}";

    return os;
}


/////////////////////////////////////////////////////////////////////////////////////////////////

/** Integration state along length of robot with Christoffel-consistent C matrix. */
struct ContinuumRobotStateChristoffel
{
    /** Number of generalized coordinates */
    unsigned n;
    /** Pose */
    Mat4r T;
    /** Pose Jacobian (6xn) */
    StrainBasisMat J;
    /** dJdot qdot / dqdot (6xn) */
    StrainBasisMat A;
    /** Mass matrix (nxn) */
    MatXr M;
    /** Christoffel-consistent Coriolis matrix, (nxn) */
    MatXr C;
    /** Generalized force (nx1) */
    VecXr Q;
    /** Potential energy */
    Real PE;
    /** Kinetic energy */
    Real KE;

    ContinuumRobotStateChristoffel() = default;

    ContinuumRobotStateChristoffel(unsigned n_)
    {
        n = n_;
        T = Mat4r::Identity();
        J = StrainBasisMat::Zero(6, n);
        A = StrainBasisMat::Zero(6, n);
        M = MatXr::Zero(n, n);
        C = MatXr::Zero(n, n);
        Q = VecXr::Zero(n);
        PE = 0;
        KE = 0;
    }

    ContinuumRobotStateChristoffel(
        unsigned n_,
        const Mat4r& T_,
        const StrainBasisMat& J_,
        const StrainBasisMat& A_,
        const MatXr& M_,
        const MatXr& C_,
        const VecXr& Q_,
        Real PE_,
        Real KE_)
        : n(n_), T(T_), J(J_), A(A_),
        M(M_), C(C_), Q(Q_), PE(PE_), KE(KE_)
    {}

    /** Overload operators for integration with RK4 */
    ContinuumRobotStateChristoffel operator+(const ContinuumRobotStateChristoffel& o) const
    {
        assert(n == o.n);
        return { n, T + o.T, J + o.J, A + o.A, M + o.M, C + o.C, Q + o.Q, PE + o.PE, KE + o.KE };
    }

    ContinuumRobotStateChristoffel operator-(const ContinuumRobotStateChristoffel& o) const
    {
        assert(n == o.n);
        return { n, T - o.T, J - o.J, A - o.A, M - o.M, C - o.C, Q - o.Q, PE - o.PE, KE - o.KE };
    }

    ContinuumRobotStateChristoffel operator*(Real s) const
    {
        return { n, s*T, s*J, s*A, s*M, s*C, s*Q, s*PE, s*KE };
    }

    ContinuumRobotStateChristoffel operator/(Real s) const
    {
        return { n, T/s, J/s, A/s, M/s, C/s, Q/s, PE/s, KE/s };
    }
};

ContinuumRobotStateChristoffel operator*(Real s, const ContinuumRobotStateChristoffel& state)
{
    return state*s;
}

std::ostream& operator<<(std::ostream& os, const ContinuumRobotStateChristoffel& state)
{
    os << "ContinuumRobotStateChristoffel {\n"
       << "  n  = " << state.n << "\n"
       << "  T  =\n" << state.T << "\n"
       << "  J  =\n" << state.J << "\n"
       << "  A  =\n" << state.A << "\n"
       << "  M  =\n" << state.M << "\n"
       << "  C  =\n" << state.C << "\n"
       << "  Q  =\n" << state.Q.transpose() << "\n"
       << "  PE = " << state.PE << "\n"
       << "  KE = " << state.KE << "\n"
       << "}";

    return os;
}