#pragma once

#include "common.hpp"

#include <vector>
#include <functional>

namespace Math
{

/** SO(3) and SE(3) calculations
 * 
 * NOTE: All calculations assume (linear, angular) ordering in se(3) =~ R^6 vectors.
 */


/** Computes the hat map from R^3 to so(3).
 * @param w - member of R^3 identified with so(3)
 * @returns the 3x3 matrix in so(3) corresponding to w
 */
Mat3r so3_hat(const Vec3r& w)
{
    Mat3r hat = Mat3r::Zero();
    hat(0,1) = -w[2];
    hat(0,2) = w[1];
    hat(1,0) = w[2];
    hat(1,2) = -w[0];
    hat(2,0) = -w[1];
    hat(2,1) = w[0];
    return hat;
}

/** Computes the hat map from R^6 to se(3).
 * @param xi - member of R^6 identified with se(3)
 * @returns the 4x4 matrix in se(3) corresponding to xi
 */
Mat4r se3_hat(const Vec6r& xi)
{
    Mat4r hat = Mat4r::Zero();
    hat.block<3,3>(0,0) = so3_hat(xi.tail<3>());
    hat.block<3,1>(0,3) = xi.head<3>();
    return hat;
}

/** Computes the adjoint of se(3).
 * @param xi - the argument to ad
 * @returns the 6x6 matrix representation of the adjoint
 */
Mat6r se3_ad(const Vec6r& xi)
{
    Mat3r vhat = so3_hat(xi.head<3>());
    Mat3r uhat = so3_hat(xi.tail<3>());
    Mat6r ad;
    ad.block<3,3>(0,0) = uhat;
    ad.block<3,3>(0,3) = vhat;
    ad.block<3,3>(3,0) = Mat3r::Zero();
    ad.block<3,3>(3,3) = uhat;
    return ad;
}

/** Computes the \bar{ad}, which is the unique matrix satisfying
 * 
 * -ad_Y^T X = adbar_X Y
 * 
 * @param w - the argument to adbar
 * @returns the 6x6 matrix corresponding to the linear operator adbar
 */
Mat6r se3_adbar(const Vec6r& w)
{
    Mat3r fhat = so3_hat(w.head<3>());
    Mat3r mhat = so3_hat(w.tail<3>());
    Mat6r adbar;
    adbar.block<3,3>(0,0) = Mat3r::Zero();
    adbar.block<3,3>(0,3) = -fhat;
    adbar.block<3,3>(3,0) = -fhat;
    adbar.block<3,3>(3,3) = -mhat;
    return adbar;
}

/** Integration */

/** General RK4 integration
 * @param state0 - the initial state at t[0]
 * @param t - a vector of integration points (t[0] corresponds to the time at state0)
 * @param params - additional constant parameters for the integration
 * @param ode_func - a function that, given (time, state, params) returns the time derivative of the state
 * @returns the states at each of the integration points
 */
template<typename StateType, typename ParamType>
static std::vector<StateType> RK4(const StateType& state0, const std::vector<Real>& t, const ParamType& params, std::function<StateType (Real, const StateType&, const ParamType&)> ode_func)
{
    std::vector<StateType> states(t.size());
    states[0] = state0;

    // run through each integration point
    for (unsigned i = 1; i < t.size(); i++)
    {
        const Real h = t[i] - t[i-1];
        
        // RK4 algorithm
        const StateType k1 = ode_func(t[i], states[i-1], params);
        const StateType k2 = ode_func(t[i] + 0.5*h, states[i-1] + 0.5*h*k1, params);
        const StateType k3 = ode_func(t[i] + 0.5*h, states[i-1] + 0.5*h*k2, params);
        const StateType k4 = ode_func(t[i] + h, states[i-1] + h*k3, params);

        states[i] = states[i-1] + h*(k1 + 2*k2 + 2*k3 + k4) / 6.0;
    }

    return states;
}

} // namespace Math