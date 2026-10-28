#pragma once

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <memory>
#include <cassert>


/** Universal typedefs used by the simulation */
using Real = double;

template <int N>
using VecNr = Eigen::Vector<Real, N>;

using Vec2r = Eigen::Vector<Real, 2>;
using Vec3r = Eigen::Vector<Real, 3>;
using Vec4r = Eigen::Vector<Real, 4>;
using Vec6r = Eigen::Vector<Real, 6>;
using VecXr = Eigen::Vector<Real, -1>;

using Vec3i = std::array<int, 3>;
using Vec4i = std::array<int, 4>;

using Vec1u = std::array<unsigned, 1>;
using Vec2u = std::array<unsigned, 2>;
using Vec3u = std::array<unsigned, 3>;
using Vec4u = std::array<unsigned, 4>;
using Vec5u = std::array<unsigned, 5>;

using Mat2r = Eigen::Matrix<Real, 2, 2>;
using Mat3r = Eigen::Matrix<Real, 3, 3>;
using Mat4r = Eigen::Matrix<Real, 4, 4>;
using Mat6r = Eigen::Matrix<Real, 6, 6>;
using MatXr = Eigen::Matrix<Real, -1, -1>;

using StrainBasisMat = Eigen::Matrix<Real, 6, -1>;
using StrainBasisFunc = std::function<StrainBasisMat (Real)>;

using Quaternion = Eigen::Quaternion<Real>;