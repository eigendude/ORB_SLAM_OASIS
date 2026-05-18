/**
* This file is part of ORB-SLAM3
*
* Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
* Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
*
* ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
* License as published by the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
* the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with ORB-SLAM3.
* If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NUMERIC_CHECKS_H
#define NUMERIC_CHECKS_H

#include <cmath>
#include <iostream>
#include <type_traits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "sophus/se3.hpp"

namespace ORB_SLAM3
{
namespace NumericChecks
{

template<typename Scalar>
inline typename std::enable_if<std::is_arithmetic<Scalar>::value, bool>::type IsFinite(const Scalar value)
{
    return std::isfinite(value);
}

template<typename Derived>
inline bool IsFinite(const Eigen::DenseBase<Derived> &value)
{
    return value.allFinite();
}

template<typename Derived>
inline bool IsFiniteQuaternion(const Eigen::QuaternionBase<Derived> &quaternion)
{
    return quaternion.coeffs().allFinite();
}

template<class Scalar, int Options>
inline bool IsFiniteSE3(const Sophus::SE3<Scalar, Options> &pose)
{
    return IsFiniteQuaternion(pose.unit_quaternion()) &&
           IsFinite(pose.translation()) &&
           IsFinite(pose.matrix());
}

template<typename RotationDerived, typename TranslationDerived>
inline bool IsFiniteSE3Parts(const Eigen::QuaternionBase<RotationDerived> &rotation,
                             const Eigen::MatrixBase<TranslationDerived> &translation)
{
    return IsFiniteQuaternion(rotation) && IsFinite(translation);
}

inline void LogRejectedInvalidFrameState(const char *description,
                                         const char *stage,
                                         const unsigned long frameId,
                                         const double timestamp)
{
    std::cerr << "ORB_SLAM_OASIS rejected invalid " << description << " in " << stage
              << " for frame " << frameId
              << " at t=" << timestamp << std::endl;
}

template<typename TranslationDerived>
inline void LogRejectedInvalidFrameState(const char *description,
                                         const char *stage,
                                         const unsigned long frameId,
                                         const double timestamp,
                                         const Eigen::MatrixBase<TranslationDerived> &translation)
{
    std::cerr << "ORB_SLAM_OASIS rejected invalid " << description << " in " << stage
              << " for frame " << frameId
              << " at t=" << timestamp
              << " tx=" << translation(0)
              << " ty=" << translation(1)
              << " tz=" << translation(2) << std::endl;
}

template<typename OmegaDerived>
inline void LogRejectedNonFiniteSim3RotationUpdate(const char *stage,
                                                   const unsigned long keyframe1Id,
                                                   const unsigned long keyframe2Id,
                                                   const Eigen::MatrixBase<OmegaDerived> &omega)
{
    std::cerr << "ORB_SLAM_OASIS rejected non-finite Sim3 rotation update in " << stage
              << " between keyframes " << keyframe1Id
              << " and " << keyframe2Id
              << " omega=[" << omega(0) << ", " << omega(1) << ", " << omega(2) << "]"
              << std::endl;
}

} // namespace NumericChecks
} // namespace ORB_SLAM3

#endif // NUMERIC_CHECKS_H
