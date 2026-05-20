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


#include "LocalMapping.h"

#include "Converter.h"
#include "GeometricTools.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Optimizer.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <sstream>

namespace ORB_SLAM3
{

namespace
{
constexpr double kVibaBiasJumpDiagnosticSeconds = 3.0;
constexpr long kVibaBiasJumpDiagnosticKeyframes = 20;

struct InertialUpdateSnapshot
{
  bool valid = false;
  unsigned long keyframe_id = 0;
  double timestamp = 0.0;
  double scale = 1.0;
  Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, -1.0);
  Eigen::Vector3f velocity = Eigen::Vector3f::Zero();
  Eigen::Vector3f gyro_bias = Eigen::Vector3f::Zero();
  Eigen::Vector3f accel_bias = Eigen::Vector3f::Zero();
  Sophus::SE3f pose;
  size_t map_points = 0;
};

struct AlignmentReprojectionStats
{
  int projected = 0;
  int matches = 0;
  int inliers = 0;
  double median_error = -1.0;
  double max_error = -1.0;
};

struct AlignmentTraceSnapshot
{
  bool valid = false;
  unsigned long frame_id = 0;
  unsigned long keyframe_id = 0;
  unsigned long map_id = 0;
  unsigned long ref_kf_id = 0;
  unsigned long first_kf_id = 0;
  unsigned long last_kf_id = 0;
  double scale = 1.0;
  std::string gravity_world_frame = "unknown";
  Eigen::Vector3d gravity_world = Eigen::Vector3d(0.0, 0.0, -1.0);
  Eigen::Vector3d gravity_camera = Eigen::Vector3d::Zero();
  Eigen::Vector3d gravity_body = Eigen::Vector3d::Zero();
  Eigen::Vector3d gravity_body_tcb_candidate = Eigen::Vector3d::Zero();
  Eigen::Vector3d specific_force_camera = Eigen::Vector3d::Zero();
  Eigen::Vector3d specific_force_body = Eigen::Vector3d::Zero();
  Eigen::Matrix3d Rbc = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Rcb = Eigen::Matrix3d::Identity();
  Sophus::SE3f Tcw;
  Sophus::SE3f ref_Twc;
  Sophus::SE3f first_Twc;
  Sophus::SE3f last_Twc;
  std::vector<Eigen::Vector3f> map_points;
  AlignmentReprojectionStats reprojection;
};

std::string FormatVector3Compact(const Eigen::Vector3f& v)
{
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(3) << "(" << v.x() << "," << v.y() << "," << v.z()
         << ")";
  return stream.str();
}

std::string FormatVector3Compact(const Eigen::Vector3d& v)
{
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(3) << "(" << v.x() << "," << v.y() << "," << v.z()
         << ")";
  return stream.str();
}

std::string FormatRpyDegCompact(const Eigen::Matrix3f& R)
{
  const Eigen::Vector3f rpy = R.eulerAngles(0, 1, 2) * 180.0f / static_cast<float>(M_PI);
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(1) << "(" << rpy.x() << "," << rpy.y() << "," << rpy.z()
         << ")";
  return stream.str();
}

std::string FormatMapPointSample(const std::vector<Eigen::Vector3f>& points, size_t index)
{
  if (index >= points.size())
    return "na";

  return FormatVector3Compact(points[index]);
}

const char* SignedAxisName(const Eigen::Vector3d& v)
{
  int axis = 0;
  double magnitude = std::abs(v.x());
  if (std::abs(v.y()) > magnitude)
  {
    axis = 1;
    magnitude = std::abs(v.y());
  }
  if (std::abs(v.z()) > magnitude)
    axis = 2;

  if (axis == 0)
    return v.x() >= 0.0 ? "+X" : "-X";
  if (axis == 1)
    return v.y() >= 0.0 ? "+Y" : "-Y";
  return v.z() >= 0.0 ? "+Z" : "-Z";
}

double UnitDot(const Eigen::Vector3d& a, const Eigen::Vector3d& b)
{
  if (a.norm() < 1e-9 || b.norm() < 1e-9)
    return 0.0;

  return a.normalized().dot(b.normalized());
}

double RotationAngleRad(const Eigen::Matrix3f& Rcw1, const Eigen::Matrix3f& Rcw2)
{
  const Eigen::Matrix3f dR = Rcw1 * Rcw2.transpose();
  const double cos_angle =
      std::max(-1.0, std::min(1.0, 0.5 * (static_cast<double>(dR.trace()) - 1.0)));
  return std::acos(cos_angle);
}

double GravityDeltaDeg(const Eigen::Vector3d& g1, const Eigen::Vector3d& g2)
{
  if (g1.norm() < 1e-9 || g2.norm() < 1e-9)
    return -1.0;

  const double cos_angle = std::max(-1.0, std::min(1.0, g1.normalized().dot(g2.normalized())));
  return std::acos(cos_angle) * 180.0 / M_PI;
}

AlignmentReprojectionStats ComputeAlignmentReprojectionStats(Frame& frame)
{
  AlignmentReprojectionStats stats;
  if (!frame.mpCamera || !frame.HasPose())
    return stats;

  std::vector<double> errors;
  const Sophus::SE3f Tcw = frame.GetPose();
  const int nLeft = frame.Nleft != -1 ? std::min(frame.Nleft, frame.N) : frame.N;
  for (int i = 0; i < nLeft; ++i)
  {
    if (i >= static_cast<int>(frame.mvpMapPoints.size()) ||
        i >= static_cast<int>(frame.mvKeysUn.size()))
      continue;

    MapPoint* pMP = frame.mvpMapPoints[i];
    if (!pMP || pMP->isBad())
      continue;

    ++stats.matches;
    if (i < static_cast<int>(frame.mvbOutlier.size()) && !frame.mvbOutlier[i])
      ++stats.inliers;

    const Eigen::Vector3f Pc = Tcw * pMP->GetWorldPos();
    if (Pc(2) <= 0.0f)
      continue;

    const Eigen::Vector2f uv = frame.mpCamera->project(Pc);
    if (!std::isfinite(uv(0)) || !std::isfinite(uv(1)))
      continue;

    ++stats.projected;
    const cv::KeyPoint& keypoint = frame.mvKeysUn[i];
    const double dx = static_cast<double>(keypoint.pt.x) - static_cast<double>(uv(0));
    const double dy = static_cast<double>(keypoint.pt.y) - static_cast<double>(uv(1));
    const double error = std::sqrt(dx * dx + dy * dy);
    errors.push_back(error);
    stats.max_error = std::max(stats.max_error, error);
  }

  if (!errors.empty())
  {
    std::sort(errors.begin(), errors.end());
    stats.median_error = errors[errors.size() / 2];
  }

  return stats;
}

AlignmentTraceSnapshot CaptureAlignmentTraceSnapshot(Tracking* tracker,
                                                     KeyFrame* current_keyframe,
                                                     double scale,
                                                     const Eigen::Vector3d& gravity_world,
                                                     const std::string& gravity_world_frame,
                                                     const std::vector<KeyFrame*>& keyframes)
{
  AlignmentTraceSnapshot snapshot;
  if (!tracker || !current_keyframe)
    return snapshot;

  Map* map = current_keyframe->GetMap();
  snapshot.valid = true;
  snapshot.frame_id = tracker->mCurrentFrame.mnId;
  snapshot.keyframe_id = current_keyframe->mnId;
  snapshot.map_id = map ? map->GetId() : 0;
  snapshot.scale = scale;
  snapshot.gravity_world = gravity_world;
  snapshot.gravity_world_frame = gravity_world_frame;
  snapshot.Tcw = tracker->mCurrentFrame.GetPose();
  snapshot.gravity_camera = snapshot.Tcw.rotationMatrix().cast<double>() * snapshot.gravity_world;
  snapshot.Rbc = tracker->mCurrentFrame.mImuCalib.mTbc.rotationMatrix().cast<double>();
  snapshot.Rcb = tracker->mCurrentFrame.mImuCalib.mTcb.rotationMatrix().cast<double>();
  snapshot.gravity_body = snapshot.Rbc * snapshot.gravity_camera;
  snapshot.gravity_body_tcb_candidate = snapshot.Rcb * snapshot.gravity_camera;
  snapshot.specific_force_camera = -snapshot.gravity_camera;
  snapshot.specific_force_body = -snapshot.gravity_body;
  snapshot.reprojection = ComputeAlignmentReprojectionStats(tracker->mCurrentFrame);

  KeyFrame* ref_keyframe = tracker->mCurrentFrame.mpReferenceKF
                               ? tracker->mCurrentFrame.mpReferenceKF
                               : current_keyframe;
  KeyFrame* first_keyframe = keyframes.empty() ? current_keyframe : keyframes.front();
  KeyFrame* last_keyframe = keyframes.empty() ? current_keyframe : keyframes.back();
  snapshot.ref_kf_id = ref_keyframe ? ref_keyframe->mnId : 0;
  snapshot.first_kf_id = first_keyframe ? first_keyframe->mnId : 0;
  snapshot.last_kf_id = last_keyframe ? last_keyframe->mnId : 0;
  snapshot.ref_Twc = ref_keyframe ? ref_keyframe->GetPoseInverse() : Sophus::SE3f();
  snapshot.first_Twc = first_keyframe ? first_keyframe->GetPoseInverse() : Sophus::SE3f();
  snapshot.last_Twc = last_keyframe ? last_keyframe->GetPoseInverse() : Sophus::SE3f();

  for (MapPoint* map_point : tracker->mCurrentFrame.mvpMapPoints)
  {
    if (!map_point || map_point->isBad())
      continue;

    snapshot.map_points.push_back(map_point->GetWorldPos());
    if (snapshot.map_points.size() == 3)
      break;
  }

  return snapshot;
}

void LogAlignmentTrace(const std::string& event_name,
                       const AlignmentTraceSnapshot& before,
                       const AlignmentTraceSnapshot& after,
                       const Eigen::Matrix3d& Rwg,
                       const float scale)
{
  if (!before.valid || !after.valid)
    return;

  const Sophus::SE3f beforeTwc = before.Tcw.inverse();
  const Sophus::SE3f afterTwc = after.Tcw.inverse();
  const Eigen::Matrix3f dRcw = after.Tcw.rotationMatrix() * before.Tcw.rotationMatrix().transpose();
  const Eigen::Vector3f dTwc = afterTwc.translation() - beforeTwc.translation();
  const double dRdeg =
      RotationAngleRad(after.Tcw.rotationMatrix(), before.Tcw.rotationMatrix()) * 180.0 / M_PI;
  const double dTwcDeg =
      RotationAngleRad(afterTwc.rotationMatrix(), beforeTwc.rotationMatrix()) * 180.0 / M_PI;
  const double dscale = after.scale - before.scale;
  const Eigen::Matrix3f Rwg_f = Rwg.cast<float>();
  const Eigen::Matrix3f Rgw_f = Rwg.transpose().cast<float>();
  const Eigen::Vector3d g_old_world = Rwg * Eigen::Vector3d(0.0, 0.0, -1.0);
  const Eigen::Vector3d wrong_epoch_gravity_camera =
      after.Tcw.rotationMatrix().cast<double>() * g_old_world;
  const Eigen::Vector3d wrong_epoch_gravity_body = after.Rbc * wrong_epoch_gravity_camera;
  const bool body_ok = UnitDot(after.gravity_body, -Eigen::Vector3d::UnitZ()) > 0.9;
  const bool camera_ok = UnitDot(after.gravity_camera, Eigen::Vector3d::UnitY()) > 0.9;

  std::cout << "MI align_formula:" << " function=LocalMapping::InitializeIMU"
            << " event=" << event_name << " input_gravity=g_old_world=mRwg*(0,0,-1)_gravity_accel"
            << " input_frame=pre_alignment_visual_world"
            << " target_gravity=g_new_world=(0,0,-1)_gravity_accel"
            << " target_frame=post_alignment_gravity_world"
            << " Rwg_rpy=" << FormatRpyDegCompact(Rwg_f)
            << " Rgw_rpy=" << FormatRpyDegCompact(Rgw_f)
            << " Rwg_angle=" << RotationAngleRad(Rwg_f, Eigen::Matrix3f::Identity()) * 180.0 / M_PI
            << " applied=Twg=Rwg^T,ApplyScaledRotation(Twg,s,true)"
            << " Tbc_applied_for_body_check=1 Tcb_applied=0"
            << " vector_semantics=gravity_accel_not_specific_force s=" << scale << std::endl;

  std::cout << "MI align_candidate:" << " event=" << event_name << " s_before=" << before.scale
            << " s_after=" << after.scale << " dscale=" << dscale
            << " dR_rpy=" << FormatRpyDegCompact(dRcw) << " dR_angle=" << dRdeg
            << " dt=" << FormatVector3Compact(dTwc) << " gravity_before_"
            << before.gravity_world_frame << "=" << FormatVector3Compact(before.gravity_world)
            << " gravity_after_" << after.gravity_world_frame << "="
            << FormatVector3Compact(after.gravity_world)
            << " gravity_before_camera=" << FormatVector3Compact(before.gravity_camera)
            << " gravity_after_camera=" << FormatVector3Compact(after.gravity_camera)
            << " gravity_before_body=" << FormatVector3Compact(before.gravity_body)
            << " gravity_after_body=" << FormatVector3Compact(after.gravity_body) << std::endl;

  std::cout << "MI align_invariant:" << " event=" << event_name
            << " gravity_camera_axis=" << SignedAxisName(after.gravity_camera)
            << " gravity_body_axis=" << SignedAxisName(after.gravity_body)
            << " specific_force_camera_axis=" << SignedAxisName(after.specific_force_camera)
            << " specific_force_body_axis=" << SignedAxisName(after.specific_force_body)
            << " gravity_camera_dot_y=" << UnitDot(after.gravity_camera, Eigen::Vector3d::UnitY())
            << " gravity_body_dot_neg_z=" << UnitDot(after.gravity_body, -Eigen::Vector3d::UnitZ())
            << " specific_force_camera_dot_neg_y="
            << UnitDot(after.specific_force_camera, -Eigen::Vector3d::UnitY())
            << " specific_force_body_dot_z="
            << UnitDot(after.specific_force_body, Eigen::Vector3d::UnitZ())
            << " ok=" << (body_ok ? 1 : 0) << std::endl;
  if (!body_ok)
  {
    std::cout << "MI align_invariant_warn:" << " event=" << event_name
              << " gravity_body_axis=" << SignedAxisName(after.gravity_body)
              << " expected=-Z possible_epoch_or_transform_mismatch=1"
              << " gravity_body=" << FormatVector3Compact(after.gravity_body)
              << " gravity_camera=" << FormatVector3Compact(after.gravity_camera)
              << " camera_ok=" << (camera_ok ? 1 : 0) << std::endl;
  }

  std::cout << "MI align_candidates:" << " event=" << event_name
            << " current_g_body_Tbc=" << FormatVector3Compact(after.gravity_body) << "/"
            << SignedAxisName(after.gravity_body)
            << " old_epoch_g_body_Tbc=" << FormatVector3Compact(wrong_epoch_gravity_body) << "/"
            << SignedAxisName(wrong_epoch_gravity_body) << " current_g_body_Tcb_wrong="
            << FormatVector3Compact(after.gravity_body_tcb_candidate) << "/"
            << SignedAxisName(after.gravity_body_tcb_candidate)
            << " specific_force_as_gravity_body_wrong="
            << FormatVector3Compact(after.specific_force_body) << "/"
            << SignedAxisName(after.specific_force_body)
            << " Rwg_maps_new_minus_z_to_old_g=1 Rgw_applied_to_map=1" << std::endl;

  std::cout << "MI align_pose:" << " event=" << event_name
            << " Tcw_before_rpy=" << FormatRpyDegCompact(before.Tcw.rotationMatrix())
            << " Tcw_after_rpy=" << FormatRpyDegCompact(after.Tcw.rotationMatrix())
            << " Twc_before_rpy=" << FormatRpyDegCompact(beforeTwc.rotationMatrix())
            << " Twc_after_rpy=" << FormatRpyDegCompact(afterTwc.rotationMatrix())
            << " dTcw_angle=" << dRdeg << " dTwc_angle=" << dTwcDeg
            << " dTwc_xyz=" << FormatVector3Compact(dTwc) << " frame_id=" << after.frame_id
            << " keyframe_id=" << after.keyframe_id << " map_id=" << after.map_id << std::endl;

  std::cout << "MI align_map_sample:" << " event=" << event_name << " ref_kf=" << after.ref_kf_id
            << " first_kf=" << after.first_kf_id << " last_kf=" << after.last_kf_id
            << " ref_Twc_before_rpy=" << FormatRpyDegCompact(before.ref_Twc.rotationMatrix())
            << " ref_Twc_after_rpy=" << FormatRpyDegCompact(after.ref_Twc.rotationMatrix())
            << " first_Twc_before_rpy=" << FormatRpyDegCompact(before.first_Twc.rotationMatrix())
            << " first_Twc_after_rpy=" << FormatRpyDegCompact(after.first_Twc.rotationMatrix())
            << " last_Twc_before_rpy=" << FormatRpyDegCompact(before.last_Twc.rotationMatrix())
            << " last_Twc_after_rpy=" << FormatRpyDegCompact(after.last_Twc.rotationMatrix())
            << " mp0_before=" << FormatMapPointSample(before.map_points, 0)
            << " mp0_after=" << FormatMapPointSample(after.map_points, 0)
            << " mp1_before=" << FormatMapPointSample(before.map_points, 1)
            << " mp1_after=" << FormatMapPointSample(after.map_points, 1)
            << " mp2_before=" << FormatMapPointSample(before.map_points, 2)
            << " mp2_after=" << FormatMapPointSample(after.map_points, 2) << std::endl;

  std::cout << "MI align_reproject:" << " event=" << event_name
            << " before_proj=" << before.reprojection.projected
            << " before_matches=" << before.reprojection.matches
            << " before_inliers=" << before.reprojection.inliers
            << " after_proj=" << after.reprojection.projected
            << " after_matches=" << after.reprojection.matches
            << " after_inliers=" << after.reprojection.inliers
            << " before_median_err=" << before.reprojection.median_error
            << " after_median_err=" << after.reprojection.median_error
            << " before_max_err=" << before.reprojection.max_error
            << " after_max_err=" << after.reprojection.max_error << std::endl;
}

InertialUpdateSnapshot CaptureInertialUpdateSnapshot(LocalMapping* local_mapper, KeyFrame* keyframe)
{
  InertialUpdateSnapshot snapshot;
  snapshot.valid = keyframe != NULL;
  snapshot.scale = local_mapper->mScale;
  snapshot.gravity = local_mapper->mRwg * Eigen::Vector3d(0.0, 0.0, -1.0);

  if (keyframe != NULL)
  {
    snapshot.keyframe_id = keyframe->mnId;
    snapshot.timestamp = keyframe->mTimeStamp;
    if (keyframe->GetMap())
      snapshot.map_points = keyframe->GetMap()->MapPointsInMap();
    snapshot.velocity = keyframe->GetVelocity();
    snapshot.gyro_bias = keyframe->GetGyroBias();
    snapshot.accel_bias = keyframe->GetAccBias();
    snapshot.pose = keyframe->GetPose();
  }

  return snapshot;
}

void LogInertialUpdateEvent(const std::string& event_name,
                            const InertialUpdateSnapshot& before,
                            const InertialUpdateSnapshot& after,
                            const bool map_points_transformed)
{
  if (!before.valid || !after.valid)
    return;

  const Sophus::SE3f beforeTwc = before.pose.inverse();
  const Sophus::SE3f afterTwc = after.pose.inverse();
  const Eigen::Vector3f dxyz = afterTwc.translation() - beforeTwc.translation();
  const double dr_deg =
      RotationAngleRad(beforeTwc.rotationMatrix(), afterTwc.rotationMatrix()) * 180.0 / M_PI;

  if (event_name == "init_accepted" || event_name == "viba_1")
  {
    std::cout << "MI inertial_event:" << " ev=" << event_name << " s=" << before.scale << "->"
              << after.scale << " dg=" << GravityDeltaDeg(before.gravity, after.gravity)
              << " dr=" << dr_deg << " dpos=" << dxyz.norm() << " pts=" << before.map_points << "->"
              << after.map_points << " mp_xform=" << map_points_transformed
              << " ba=" << after.accel_bias.norm() << " bg=" << after.gyro_bias.norm() << std::endl;
  }
}

void LogLocalInertialBA(const InertialUpdateSnapshot& before,
                        const InertialUpdateSnapshot& after,
                        const Optimizer::LocalInertialBADiagnostic& diagnostic,
                        const double since_viba_sec,
                        const long since_viba_kfs)
{
  if (!before.valid || !after.valid)
    return;

  const Sophus::SE3f beforeTwc = before.pose.inverse();
  const Sophus::SE3f afterTwc = after.pose.inverse();
  const Eigen::Vector3f dxyz = afterTwc.translation() - beforeTwc.translation();
  const double dr_deg =
      RotationAngleRad(beforeTwc.rotationMatrix(), afterTwc.rotationMatrix()) * 180.0 / M_PI;
  const Eigen::Vector3f dbg = after.gyro_bias - before.gyro_bias;
  const Eigen::Vector3f dba = after.accel_bias - before.accel_bias;
  const bool map_drop = before.map_points > 0 && after.map_points * 10 < before.map_points * 9;
  static double last_log_time = -1e9;
  const bool enough_time = after.timestamp - last_log_time >= 1.0;
  const bool significant = dba.norm() > 0.03f || after.accel_bias.norm() > 0.12f || dr_deg > 1.0 ||
                           dxyz.norm() > 0.05f || map_drop;

  if (enough_time || significant)
  {
    last_log_time = after.timestamp;
    std::cout << "MI iba:" << " s=" << after.scale << " kfs=" << diagnostic.optimized_keyframes
              << "/" << diagnostic.fixed_keyframes << " edges=" << diagnostic.visual_edges << "/"
              << diagnostic.inertial_edges << "/" << diagnostic.bias_random_walk_edges
              << " pts=" << before.map_points << "->" << after.map_points
              << " ba=" << after.accel_bias.norm() << " dba=" << dba.norm()
              << " bg=" << after.gyro_bias.norm() << " dbg=" << dbg.norm()
              << " chi_v=" << diagnostic.visual_chi2_before << "->" << diagnostic.visual_chi2_after
              << " chi_i=" << diagnostic.inertial_chi2_before << "->"
              << diagnostic.inertial_chi2_after << " dr=" << dr_deg << " dpos=" << dxyz.norm()
              << " vel=" << after.velocity.norm() << " bias_src=current_kf_optimizer" << std::endl;
  }

  const bool post_viba_diag_window =
      since_viba_sec >= 0.0 && (since_viba_sec <= kVibaBiasJumpDiagnosticSeconds ||
                                since_viba_kfs <= kVibaBiasJumpDiagnosticKeyframes);
  if (post_viba_diag_window && dba.norm() > 0.05f)
  {
    std::cout << "MI ba_jump:" << " since_viba=" << since_viba_sec << "s/" << since_viba_kfs << "kf"
              << " ba_before=" << FormatVector3Compact(before.accel_bias)
              << " ba_after=" << FormatVector3Compact(after.accel_bias) << " dba=" << dba.norm()
              << " inertial_chi2=" << diagnostic.inertial_chi2_before << "->"
              << diagnostic.inertial_chi2_after << " visual_chi2=" << diagnostic.visual_chi2_before
              << "->" << diagnostic.visual_chi2_after << " kfs=" << diagnostic.optimized_keyframes
              << "/" << diagnostic.fixed_keyframes << " edges=" << diagnostic.visual_edges << "/"
              << diagnostic.inertial_edges << "/" << diagnostic.bias_random_walk_edges
              << " curve=na" << std::endl;
  }
}

} // namespace

LocalMapping::LocalMapping(System* pSys, Atlas *pAtlas, const float bMonocular, bool bInertial, const string &_strSeqName):
    mpSystem(pSys), mbMonocular(bMonocular), mbInertial(bInertial), mbResetRequested(false), mbResetRequestedActiveMap(false), mbFinishRequested(false), mbFinished(true), mpAtlas(pAtlas), bInitializing(false),
    mbAbortBA(false), mbStopped(false), mbStopRequested(false), mbNotStop(false), mbAcceptKeyFrames(true),
    mIdxInit(0), mScale(1.0), mInitSect(0), mbNotBA1(true), mbNotBA2(true), mIdxIteration(0), infoInertial(Eigen::MatrixXd::Zero(9,9))
{
    mnMatchesInliers = 0;

    mbBadImu = false;
    mLastViba1Time = -1.0;
    mLastViba1KfId = 0;
    mInitialImuKeyframeId = 0;
    mInertialInitAttemptId = 0;
    mLastInitLifecycleCommitted = false;
    mLastAlignmentEvent = "none";
    mLastAlignmentRotationDeg = 0.0;
    mLastAlignmentScale = 1.0;

    mTinit = 0.f;

    mNumLM = 0;
    mNumKFCulling=0;

#ifdef REGISTER_TIMES
    nLBA_exec = 0;
    nLBA_abort = 0;
#endif

}

bool LocalMapping::IsInertialInitializationProvisional() const
{
  if (!mbInertial || !mbMonocular || mbBadImu || mpAtlas == nullptr)
    return false;

  Map* pMap = mpAtlas->GetCurrentMap();
  if (pMap == nullptr || !pMap->isImuInitialized())
    return false;

  return !pMap->GetIniertialBA2() && mTinit < 10.0f;
}

bool LocalMapping::IsInertialInitializationCommitted() const
{
  if (!mbInertial || mpAtlas == nullptr)
    return false;

  Map* pMap = mpAtlas->GetCurrentMap();
  if (pMap == nullptr || !pMap->isImuInitialized() || mbBadImu)
    return false;

  return !IsInertialInitializationProvisional();
}

void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

void LocalMapping::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}

void LocalMapping::Run()
{
    mbFinished = false;

    while(1)
    {
        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(false);

        // Check if there are keyframes in the queue
        if(CheckNewKeyFrames() && !mbBadImu)
        {
#ifdef REGISTER_TIMES
            double timeLBA_ms = 0;
            double timeKFCulling_ms = 0;

            std::chrono::steady_clock::time_point time_StartProcessKF = std::chrono::steady_clock::now();
#endif
            // BoW conversion and insertion in Map
            ProcessNewKeyFrame();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndProcessKF = std::chrono::steady_clock::now();

            double timeProcessKF = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndProcessKF - time_StartProcessKF).count();
            vdKFInsert_ms.push_back(timeProcessKF);
#endif

            // Check recent MapPoints
            MapPointCulling();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCulling = std::chrono::steady_clock::now();

            double timeMPCulling = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMPCulling - time_EndProcessKF).count();
            vdMPCulling_ms.push_back(timeMPCulling);
#endif

            // Triangulate new MapPoints
            CreateNewMapPoints();

            mbAbortBA = false;

            if(!CheckNewKeyFrames())
            {
                // Find more matches in neighbor keyframes and fuse point duplications
                SearchInNeighbors();
            }

#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCreation = std::chrono::steady_clock::now();

            double timeMPCreation = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndMPCreation - time_EndMPCulling).count();
            vdMPCreation_ms.push_back(timeMPCreation);
#endif

            bool b_doneLBA = false;
            int num_FixedKF_BA = 0;
            int num_OptKF_BA = 0;
            int num_MPs_BA = 0;
            int num_edges_BA = 0;

            if(!CheckNewKeyFrames() && !stopRequested())
            {
                if(mpAtlas->KeyFramesInMap()>2)
                {

                    if(mbInertial && mpCurrentKeyFrame->GetMap()->isImuInitialized())
                    {
                        float dist = (mpCurrentKeyFrame->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->GetCameraCenter()).norm() +
                                (mpCurrentKeyFrame->mPrevKF->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->mPrevKF->GetCameraCenter()).norm();

                        if(dist>0.05)
                            mTinit += mpCurrentKeyFrame->mTimeStamp - mpCurrentKeyFrame->mPrevKF->mTimeStamp;
                        const bool motion_check_not_enough =
                            !mpCurrentKeyFrame->GetMap()->GetIniertialBA2() && (mTinit < 10.f) &&
                            (dist < 0.02);
                        if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2())
                        {
                          if (motion_check_not_enough)
                          {
                            std::cout << "MI init_rejected:" << " reason=not_enough_motion"
                                      << " mTinit=" << mTinit << " dist=" << dist
                                      << " kfs=" << mpAtlas->KeyFramesInMap()
                                      << " map=" << mpAtlas->MapPointsInMap()
                                      << " inl=" << mpTracker->GetMatchesInliers() << std::endl;
                            std::cout
                                << "MI init_lifecycle:" << " attempt=" << mInertialInitAttemptId
                                << " phase=motion_check after_align=1 map_imu="
                                << (mpCurrentKeyFrame->GetMap()->isImuInitialized() ? 1 : 0)
                                << " tracking_state=" << mpTracker->mState
                                << " bad_imu=1 motion=not_enough should_publish=0 committed=0"
                                << " reset=1 mTinit=" << mTinit << " dist=" << dist
                                << " kfs=" << mpAtlas->KeyFramesInMap()
                                << " map=" << mpAtlas->MapPointsInMap() << std::endl;
                            cout << "Not enough motion for initializing. Reseting..." << endl;
                            if (mpTracker)
                            {
                              mpTracker->RecordTrackingFailure(
                                  TrackingFailureReason::NotEnoughMotionForImuInitialization,
                                  mpCurrentKeyFrame->mTimeStamp);
                            }
                            std::unique_lock<std::mutex> lock(mMutexReset);
                            mbResetRequestedActiveMap = true;
                            mpMapToReset = mpCurrentKeyFrame->GetMap();
                            mbBadImu = true;
                          }
                        }
                        if (!motion_check_not_enough && IsInertialInitializationCommitted() &&
                            !mLastInitLifecycleCommitted)
                        {
                          mLastInitLifecycleCommitted = true;
                          std::cout << "MI init_lifecycle:" << " attempt=" << mInertialInitAttemptId
                                    << " phase=motion_check after_align=1 map_imu=1"
                                    << " tracking_state=" << mpTracker->mState
                                    << " bad_imu=0 motion=ok should_publish=1 committed=1 reset=0"
                                    << " mTinit=" << mTinit << " dist=" << dist
                                    << " kfs=" << mpAtlas->KeyFramesInMap()
                                    << " map=" << mpAtlas->MapPointsInMap() << std::endl;
                        }

                        bool bLarge = ((mpTracker->GetMatchesInliers()>75)&&mbMonocular)||((mpTracker->GetMatchesInliers()>100)&&!mbMonocular);
                        const InertialUpdateSnapshot before_inertial_ba =
                            CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame);
                        Optimizer::LocalInertialBADiagnostic iba_diagnostic;
                        const double since_viba_sec =
                            mLastViba1Time >= 0.0 ? mpCurrentKeyFrame->mTimeStamp - mLastViba1Time
                                                  : -1.0;
                        const long since_viba_kfs =
                            mLastViba1Time >= 0.0 ? static_cast<long>(mpCurrentKeyFrame->mnId) -
                                                        static_cast<long>(mLastViba1KfId)
                                                  : -1;
                        Optimizer::LocalInertialBA(
                            mpCurrentKeyFrame, &mbAbortBA, mpCurrentKeyFrame->GetMap(),
                            num_FixedKF_BA, num_OptKF_BA, num_MPs_BA, num_edges_BA, bLarge,
                            !mpCurrentKeyFrame->GetMap()->GetIniertialBA2(), &iba_diagnostic);
                        LogLocalInertialBA(before_inertial_ba,
                                           CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame),
                                           iba_diagnostic, since_viba_sec, since_viba_kfs);
                        b_doneLBA = true;
                    }
                    else
                    {
                        Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA, mpCurrentKeyFrame->GetMap(),num_FixedKF_BA,num_OptKF_BA,num_MPs_BA,num_edges_BA);
                        b_doneLBA = true;
                    }

                }
#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndLBA = std::chrono::steady_clock::now();

                if(b_doneLBA)
                {
                    timeLBA_ms = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLBA - time_EndMPCreation).count();
                    vdLBA_ms.push_back(timeLBA_ms);

                    nLBA_exec += 1;
                    if(mbAbortBA)
                    {
                        nLBA_abort += 1;
                    }
                    vnLBA_edges.push_back(num_edges_BA);
                    vnLBA_KFopt.push_back(num_OptKF_BA);
                    vnLBA_KFfixed.push_back(num_FixedKF_BA);
                    vnLBA_MPs.push_back(num_MPs_BA);
                }

#endif

                // Initialize IMU here
                if(!mpCurrentKeyFrame->GetMap()->isImuInitialized() && mbInertial)
                {
                  const InertialUpdateSnapshot before_init =
                      CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame);
                  if (mbMonocular)
                    InitializeIMU(1e2, 1e10, true);
                  else
                    InitializeIMU(1e2, 1e5, true);
                  if (mpCurrentKeyFrame->GetMap()->isImuInitialized())
                  {
                    LogInertialUpdateEvent("init_accepted", before_init,
                                           CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame),
                                           true);
                  }
                }


                // Check redundant local Keyframes
                KeyFrameCulling();

#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndKFCulling = std::chrono::steady_clock::now();

                timeKFCulling_ms = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndKFCulling - time_EndLBA).count();
                vdKFCulling_ms.push_back(timeKFCulling_ms);
#endif

                if ((mTinit<50.0f) && mbInertial)
                {
                    if(mpCurrentKeyFrame->GetMap()->isImuInitialized() && mpTracker->mState==Tracking::OK) // Enter here everytime local-mapping is called
                    {
                        if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA1()){
                            if (mTinit>5.0f)
                            {
                                cout << "start VIBA 1" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA1();
                                const InertialUpdateSnapshot before_viba =
                                    CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame);
                                if (mbMonocular)
                                    InitializeIMU(1.f, 1e5, true);
                                else
                                    InitializeIMU(1.f, 1e5, true);
                                LogInertialUpdateEvent(
                                    "viba_1", before_viba,
                                    CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame), true);
                                mLastViba1Time = mpCurrentKeyFrame->mTimeStamp;
                                mLastViba1KfId = mpCurrentKeyFrame->mnId;

                                cout << "end VIBA 1" << endl;
                            }
                        }
                        else if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2()){
                            if (mTinit>15.0f){
                                cout << "start VIBA 2" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA2();
                                const InertialUpdateSnapshot before_viba =
                                    CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame);
                                if (mbMonocular)
                                    InitializeIMU(0.f, 0.f, true);
                                else
                                    InitializeIMU(0.f, 0.f, true);
                                LogInertialUpdateEvent(
                                    "viba_2", before_viba,
                                    CaptureInertialUpdateSnapshot(this, mpCurrentKeyFrame), true);

                                cout << "end VIBA 2" << endl;
                            }
                        }

                        // scale refinement
                        if (((mpAtlas->KeyFramesInMap())<=200) &&
                                ((mTinit>25.0f && mTinit<25.5f)||
                                (mTinit>35.0f && mTinit<35.5f)||
                                (mTinit>45.0f && mTinit<45.5f)||
                                (mTinit>55.0f && mTinit<55.5f)||
                                (mTinit>65.0f && mTinit<65.5f)||
                                (mTinit>75.0f && mTinit<75.5f))){
                            if (mbMonocular)
                                ScaleRefinement();
                        }
                    }
                }
            }

#ifdef REGISTER_TIMES
            vdLBASync_ms.push_back(timeKFCulling_ms);
            vdKFCullingSync_ms.push_back(timeKFCulling_ms);
#endif

            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);

#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndLocalMap = std::chrono::steady_clock::now();

            double timeLocalMap = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(time_EndLocalMap - time_StartProcessKF).count();
            vdLMTotal_ms.push_back(timeLocalMap);
#endif
        }

        else if(Stop() && !mbBadImu)
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                usleep(3000);
            }
            if(CheckFinish())
                break;
        }

        ResetIfRequested();

        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(true);

        if(CheckFinish())
            break;

        usleep(3000);
    }

    SetFinish();
}

void LocalMapping::InsertKeyFrame(KeyFrame *pKF)
{
    std::unique_lock<std::mutex> lock(mMutexNewKFs);
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA=true;
}


bool LocalMapping::CheckNewKeyFrames()
{
    std::unique_lock<std::mutex> lock(mMutexNewKFs);
    return(!mlNewKeyFrames.empty());
}

void LocalMapping::ProcessNewKeyFrame()
{
    {
        std::unique_lock<std::mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    pMP->UpdateNormalAndDepth();
                    pMP->ComputeDistinctiveDescriptors();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }

    // Update links in the Covisibility Graph
    mpCurrentKeyFrame->UpdateConnections();

    // Insert Keyframe in Map
    mpAtlas->AddKeyFrame(mpCurrentKeyFrame);
}

void LocalMapping::EmptyQueue()
{
    while(CheckNewKeyFrames())
        ProcessNewKeyFrame();
}

void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    int borrar = mlpRecentAddedMapPoints.size();

    while(lit!=mlpRecentAddedMapPoints.end())
    {
        MapPoint* pMP = *lit;

        if(pMP->isBad())
            lit = mlpRecentAddedMapPoints.erase(lit);
        else if(pMP->GetFoundRatio()<0.25f)
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=2 && pMP->Observations()<=cnThObs)
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=3)
            lit = mlpRecentAddedMapPoints.erase(lit);
        else
        {
            lit++;
            borrar--;
        }
    }
}


void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    // For stereo inertial case
    if(mbMonocular)
        nn=30;
    vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    if (mbInertial)
    {
        KeyFrame* pKF = mpCurrentKeyFrame;
        int count=0;
        while((vpNeighKFs.size()<=nn)&&(pKF->mPrevKF)&&(count++<nn))
        {
            vector<KeyFrame*>::iterator it = std::find(vpNeighKFs.begin(), vpNeighKFs.end(), pKF->mPrevKF);
            if(it==vpNeighKFs.end())
                vpNeighKFs.push_back(pKF->mPrevKF);
            pKF = pKF->mPrevKF;
        }
    }

    float th = 0.6f;

    ORBmatcher matcher(th,false);

    Sophus::SE3<float> sophTcw1 = mpCurrentKeyFrame->GetPose();
    Eigen::Matrix<float,3,4> eigTcw1 = sophTcw1.matrix3x4();
    Eigen::Matrix<float,3,3> Rcw1 = eigTcw1.block<3,3>(0,0);
    Eigen::Matrix<float,3,3> Rwc1 = Rcw1.transpose();
    Eigen::Vector3f tcw1 = sophTcw1.translation();
    Eigen::Vector3f Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactor;
    int countStereo = 0;
    int countStereoGoodProj = 0;
    int countStereoAttempt = 0;
    int totalStereoPts = 0;
    // Search matches with epipolar restriction and triangulate
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {
        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        GeometricCamera* pCamera1 = mpCurrentKeyFrame->mpCamera, *pCamera2 = pKF2->mpCamera;

        // Check first that baseline is not too short
        Eigen::Vector3f Ow2 = pKF2->GetCameraCenter();
        Eigen::Vector3f vBaseline = Ow2-Ow1;
        const float baseline = vBaseline.norm();

        if(!mbMonocular)
        {
            if(baseline<pKF2->mb)
                continue;
        }
        else
        {
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            const float ratioBaselineDepth = baseline/medianDepthKF2;

            if(ratioBaselineDepth<0.01)
                continue;
        }

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        bool bCoarse = mbInertial && mpTracker->mState==Tracking::RECENTLY_LOST && mpCurrentKeyFrame->GetMap()->GetIniertialBA2();

        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,vMatchedIndices,false,bCoarse);

        Sophus::SE3<float> sophTcw2 = pKF2->GetPose();
        Eigen::Matrix<float,3,4> eigTcw2 = sophTcw2.matrix3x4();
        Eigen::Matrix<float,3,3> Rcw2 = eigTcw2.block<3,3>(0,0);
        Eigen::Matrix<float,3,3> Rwc2 = Rcw2.transpose();
        Eigen::Vector3f tcw2 = sophTcw2.translation();

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = (mpCurrentKeyFrame -> NLeft == -1) ? mpCurrentKeyFrame->mvKeysUn[idx1]
                                                                         : (idx1 < mpCurrentKeyFrame -> NLeft) ? mpCurrentKeyFrame -> mvKeys[idx1]
                                                                                                               : mpCurrentKeyFrame -> mvKeysRight[idx1 - mpCurrentKeyFrame -> NLeft];
            const float kp1_ur=mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = (!mpCurrentKeyFrame->mpCamera2 && kp1_ur>=0);
            const bool bRight1 = (mpCurrentKeyFrame -> NLeft == -1 || idx1 < mpCurrentKeyFrame -> NLeft) ? false
                                                                                                         : true;

            const cv::KeyPoint &kp2 = (pKF2 -> NLeft == -1) ? pKF2->mvKeysUn[idx2]
                                                            : (idx2 < pKF2 -> NLeft) ? pKF2 -> mvKeys[idx2]
                                                                                     : pKF2 -> mvKeysRight[idx2 - pKF2 -> NLeft];

            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = (!pKF2->mpCamera2 && kp2_ur>=0);
            const bool bRight2 = (pKF2 -> NLeft == -1 || idx2 < pKF2 -> NLeft) ? false
                                                                               : true;

            if(mpCurrentKeyFrame->mpCamera2 && pKF2->mpCamera2){
                if(bRight1 && bRight2){
                    sophTcw1 = mpCurrentKeyFrame->GetRightPose();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter();

                    sophTcw2 = pKF2->GetRightPose();
                    Ow2 = pKF2->GetRightCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera2;
                }
                else if(bRight1 && !bRight2){
                    sophTcw1 = mpCurrentKeyFrame->GetRightPose();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter();

                    sophTcw2 = pKF2->GetPose();
                    Ow2 = pKF2->GetCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera;
                }
                else if(!bRight1 && bRight2){
                    sophTcw1 = mpCurrentKeyFrame->GetPose();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter();

                    sophTcw2 = pKF2->GetRightPose();
                    Ow2 = pKF2->GetRightCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera2;
                }
                else{
                    sophTcw1 = mpCurrentKeyFrame->GetPose();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter();

                    sophTcw2 = pKF2->GetPose();
                    Ow2 = pKF2->GetCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera;
                }
                eigTcw1 = sophTcw1.matrix3x4();
                Rcw1 = eigTcw1.block<3,3>(0,0);
                Rwc1 = Rcw1.transpose();
                tcw1 = sophTcw1.translation();

                eigTcw2 = sophTcw2.matrix3x4();
                Rcw2 = eigTcw2.block<3,3>(0,0);
                Rwc2 = Rcw2.transpose();
                tcw2 = sophTcw2.translation();
            }

            // Check parallax between rays
            Eigen::Vector3f xn1 = pCamera1->unprojectEig(kp1.pt);
            Eigen::Vector3f xn2 = pCamera2->unprojectEig(kp2.pt);

            Eigen::Vector3f ray1 = Rwc1 * xn1;
            Eigen::Vector3f ray2 = Rwc2 * xn2;
            const float cosParallaxRays = ray1.dot(ray2)/(ray1.norm() * ray2.norm());

            float cosParallaxStereo = cosParallaxRays+1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            if(bStereo1)
                cosParallaxStereo1 = cos(2*atan2(mpCurrentKeyFrame->mb/2,mpCurrentKeyFrame->mvDepth[idx1]));
            else if(bStereo2)
                cosParallaxStereo2 = cos(2*atan2(pKF2->mb/2,pKF2->mvDepth[idx2]));

            if (bStereo1 || bStereo2) totalStereoPts++;
            
            cosParallaxStereo = min(cosParallaxStereo1,cosParallaxStereo2);

            Eigen::Vector3f x3D;

            bool goodProj = false;
            bool bPointStereo = false;
            if(cosParallaxRays<cosParallaxStereo && cosParallaxRays>0 && (bStereo1 || bStereo2 ||
                                                                          (cosParallaxRays<0.9996 && mbInertial) || (cosParallaxRays<0.9998 && !mbInertial)))
            {
                goodProj = GeometricTools::Triangulate(xn1, xn2, eigTcw1, eigTcw2, x3D);
                if(!goodProj)
                    continue;
            }
            else if(bStereo1 && cosParallaxStereo1<cosParallaxStereo2)
            {
                countStereoAttempt++;
                bPointStereo = true;
                goodProj = mpCurrentKeyFrame->UnprojectStereo(idx1, x3D);
            }
            else if(bStereo2 && cosParallaxStereo2<cosParallaxStereo1)
            {
                countStereoAttempt++;
                bPointStereo = true;
                goodProj = pKF2->UnprojectStereo(idx2, x3D);
            }
            else
            {
                continue; //No stereo and very low parallax
            }

            if(goodProj && bPointStereo)
                countStereoGoodProj++;

            if(!goodProj)
                continue;

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3D) + tcw1(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3D) + tcw2(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3D)+tcw1(0);
            const float y1 = Rcw1.row(1).dot(x3D)+tcw1(1);
            const float invz1 = 1.0/z1;

            if(!bStereo1)
            {
                cv::Point2f uv1 = pCamera1->project(cv::Point3f(x1,y1,z1));
                float errX1 = uv1.x - kp1.pt.x;
                float errY1 = uv1.y - kp1.pt.y;

                if((errX1*errX1+errY1*errY1)>5.991*sigmaSquare1)
                    continue;

            }
            else
            {
                float u1 = fx1*x1*invz1+cx1;
                float u1_r = u1 - mpCurrentKeyFrame->mbf*invz1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                if((errX1*errX1+errY1*errY1+errX1_r*errX1_r)>7.8*sigmaSquare1)
                    continue;
            }

            //Check reprojection error in second keyframe
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3D)+tcw2(0);
            const float y2 = Rcw2.row(1).dot(x3D)+tcw2(1);
            const float invz2 = 1.0/z2;
            if(!bStereo2)
            {
                cv::Point2f uv2 = pCamera2->project(cv::Point3f(x2,y2,z2));
                float errX2 = uv2.x - kp2.pt.x;
                float errY2 = uv2.y - kp2.pt.y;
                if((errX2*errX2+errY2*errY2)>5.991*sigmaSquare2)
                    continue;
            }
            else
            {
                float u2 = fx2*x2*invz2+cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf*invz2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if((errX2*errX2+errY2*errY2+errX2_r*errX2_r)>7.8*sigmaSquare2)
                    continue;
            }

            //Check scale consistency
            Eigen::Vector3f normal1 = x3D - Ow1;
            float dist1 = normal1.norm();

            Eigen::Vector3f normal2 = x3D - Ow2;
            float dist2 = normal2.norm();

            if(dist1==0 || dist2==0)
                continue;

            if(mbFarPoints && (dist1>=mThFarPoints||dist2>=mThFarPoints)) // MODIFICATION
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave]/pKF2->mvScaleFactors[kp2.octave];

            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;

            // Triangulation is succesfull
            MapPoint* pMP = new MapPoint(x3D, mpCurrentKeyFrame, mpAtlas->GetCurrentMap());
            if (bPointStereo)
                countStereo++;
            
            pMP->AddObservation(mpCurrentKeyFrame,idx1);
            pMP->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP,idx1);
            pKF2->AddMapPoint(pMP,idx2);

            pMP->ComputeDistinctiveDescriptors();

            pMP->UpdateNormalAndDepth();

            mpAtlas->AddMapPoint(pMP);
            mlpRecentAddedMapPoints.push_back(pMP);
        }
    }    
}

void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if(mbMonocular)
        nn=30;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    vector<KeyFrame*> vpTargetKFs;
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        if(pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;
    }

    // Add some covisible of covisible
    // Extend to some second neighbors if abort is not requested
    for(int i=0, imax=vpTargetKFs.size(); i<imax; i++)
    {
        const vector<KeyFrame*> vpSecondNeighKFs = vpTargetKFs[i]->GetBestCovisibilityKeyFrames(20);
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            if(pKFi2->isBad() || pKFi2->mnFuseTargetForKF==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
            pKFi2->mnFuseTargetForKF=mpCurrentKeyFrame->mnId;
        }
        if (mbAbortBA)
            break;
    }

    // Extend to temporal neighbors
    if(mbInertial)
    {
        KeyFrame* pKFi = mpCurrentKeyFrame->mPrevKF;
        while(vpTargetKFs.size()<20 && pKFi)
        {
            if(pKFi->isBad() || pKFi->mnFuseTargetForKF==mpCurrentKeyFrame->mnId)
            {
                pKFi = pKFi->mPrevKF;
                continue;
            }
            vpTargetKFs.push_back(pKFi);
            pKFi->mnFuseTargetForKF=mpCurrentKeyFrame->mnId;
            pKFi = pKFi->mPrevKF;
        }
    }

    // Search matches by projection from current KF in target KFs
    ORBmatcher matcher;
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi,vpMapPointMatches);
        if(pKFi->NLeft != -1) matcher.Fuse(pKFi,vpMapPointMatches,true);
    }


    if (mbAbortBA)
        return;

    // Search matches by projection from target KFs in current KF
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapPointMatches.size());

    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        for(vector<MapPoint*>::iterator vitMP=vpMapPointsKFi.begin(), vendMP=vpMapPointsKFi.end(); vitMP!=vendMP; vitMP++)
        {
            MapPoint* pMP = *vitMP;
            if(!pMP)
                continue;
            if(pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);
    if(mpCurrentKeyFrame->NLeft != -1) matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates,true);


    // Update points
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
        MapPoint* pMP=vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                pMP->ComputeDistinctiveDescriptors();
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    mpCurrentKeyFrame->UpdateConnections();
}

void LocalMapping::RequestStop()
{
    std::unique_lock<std::mutex> lock(mMutexStop);
    mbStopRequested = true;
    std::unique_lock<std::mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

bool LocalMapping::Stop()
{
    std::unique_lock<std::mutex> lock(mMutexStop);
    if(mbStopRequested && !mbNotStop)
    {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

bool LocalMapping::isStopped()
{
    std::unique_lock<std::mutex> lock(mMutexStop);
    return mbStopped;
}

bool LocalMapping::stopRequested()
{
    std::unique_lock<std::mutex> lock(mMutexStop);
    return mbStopRequested;
}

void LocalMapping::Release()
{
    std::unique_lock<std::mutex> lock(mMutexStop);
    std::unique_lock<std::mutex> lock2(mMutexFinish);
    if(mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

bool LocalMapping::AcceptKeyFrames()
{
    std::unique_lock<std::mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    std::unique_lock<std::mutex> lock(mMutexAccept);
    mbAcceptKeyFrames=flag;
}

bool LocalMapping::SetNotStop(bool flag)
{
    std::unique_lock<std::mutex> lock(mMutexStop);

    if(flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

void LocalMapping::KeyFrameCulling()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
    const int Nd = 21;
    mpCurrentKeyFrame->UpdateBestCovisibles();
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    float redundant_th;
    if(!mbInertial)
        redundant_th = 0.9;
    else if (mbMonocular)
        redundant_th = 0.9;
    else
        redundant_th = 0.5;

    const bool bInitImu = mpAtlas->isImuInitialized();
    int count=0;

    // Compoute last KF from optimizable window:
    unsigned int last_ID;
    if (mbInertial)
    {
        int count = 0;
        KeyFrame* aux_KF = mpCurrentKeyFrame;
        while(count<Nd && aux_KF->mPrevKF)
        {
            aux_KF = aux_KF->mPrevKF;
            count++;
        }
        last_ID = aux_KF->mnId;
    }



    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        count++;
        KeyFrame* pKF = *vit;

        if((pKF->mnId==pKF->GetMap()->GetInitKFid()) || pKF->isBad())
            continue;
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMPs=0;
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = (pKF -> NLeft == -1) ? pKF->mvKeysUn[i].octave
                                                                     : (i < pKF -> NLeft) ? pKF -> mvKeys[i].octave
                                                                                          : pKF -> mvKeysRight[i].octave;
                        const map<KeyFrame*, tuple<int,int>> observations = pMP->GetObservations();
                        int nObs=0;
                        for(map<KeyFrame*, tuple<int,int>>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            tuple<int,int> indexes = mit->second;
                            int leftIndex = get<0>(indexes), rightIndex = get<1>(indexes);
                            int scaleLeveli = -1;
                            if(pKFi -> NLeft == -1)
                                scaleLeveli = pKFi->mvKeysUn[leftIndex].octave;
                            else {
                                if (leftIndex != -1) {
                                    scaleLeveli = pKFi->mvKeys[leftIndex].octave;
                                }
                                if (rightIndex != -1) {
                                    int rightLevel = pKFi->mvKeysRight[rightIndex - pKFi->NLeft].octave;
                                    scaleLeveli = (scaleLeveli == -1 || scaleLeveli > rightLevel) ? rightLevel
                                                                                                  : scaleLeveli;
                                }
                            }

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>thObs)
                                    break;
                            }
                        }
                        if(nObs>thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }

        if(nRedundantObservations>redundant_th*nMPs)
        {
            if (mbInertial)
            {
                if (mpAtlas->KeyFramesInMap()<=Nd)
                    continue;

                if(pKF->mnId>(mpCurrentKeyFrame->mnId-2))
                    continue;

                if(pKF->mPrevKF && pKF->mNextKF)
                {
                    const float t = pKF->mNextKF->mTimeStamp-pKF->mPrevKF->mTimeStamp;

                    if((bInitImu && (pKF->mnId<last_ID) && t<3.) || (t<0.5))
                    {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    }
                    else if(!mpCurrentKeyFrame->GetMap()->GetIniertialBA2() && ((pKF->GetImuPosition()-pKF->mPrevKF->GetImuPosition()).norm()<0.02) && (t<3))
                    {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    }
                }
            }
            else
            {
                pKF->SetBadFlag();
            }
        }
        if((count > 20 && mbAbortBA) || count>100)
        {
            break;
        }
    }
}

void LocalMapping::RequestReset()
{
    {
        std::unique_lock<std::mutex> lock(mMutexReset);
        cout << "LM: Map reset recieved" << endl;
        mbResetRequested = true;
    }
    cout << "LM: Map reset, waiting..." << endl;

    while(1)
    {
        {
            std::unique_lock<std::mutex> lock2(mMutexReset);
            if(!mbResetRequested)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Map reset, Done!!!" << endl;
}

void LocalMapping::RequestResetActiveMap(Map* pMap)
{
    {
        std::unique_lock<std::mutex> lock(mMutexReset);
        cout << "LM: Active map reset recieved" << endl;
        mbResetRequestedActiveMap = true;
        mpMapToReset = pMap;
    }
    cout << "LM: Active map reset, waiting..." << endl;

    while(1)
    {
        {
            std::unique_lock<std::mutex> lock2(mMutexReset);
            if(!mbResetRequestedActiveMap)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Active map reset, Done!!!" << endl;
}

void LocalMapping::ResetIfRequested()
{
    bool executed_reset = false;
    {
        std::unique_lock<std::mutex> lock(mMutexReset);
        if(mbResetRequested)
        {
            executed_reset = true;

            cout << "LM: Reseting Atlas in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();
            mbResetRequested = false;
            mbResetRequestedActiveMap = false;

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu=false;
            mLastInitLifecycleCommitted = false;

            mIdxInit=0;

            cout << "LM: End reseting Local Mapping..." << endl;
        }

        if(mbResetRequestedActiveMap) {
            executed_reset = true;
            cout << "LM: Reseting current map in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu=false;
            mLastInitLifecycleCommitted = false;

            mbResetRequested = false;
            mbResetRequestedActiveMap = false;
            cout << "LM: End reseting Local Mapping..." << endl;
        }
    }
    if(executed_reset)
        cout << "LM: Reset free the mutex" << endl;

}

void LocalMapping::RequestFinish()
{
    std::unique_lock<std::mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

bool LocalMapping::CheckFinish()
{
    std::unique_lock<std::mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

void LocalMapping::SetFinish()
{
    std::unique_lock<std::mutex> lock(mMutexFinish);
    mbFinished = true;    
    std::unique_lock<std::mutex> lock2(mMutexStop);
    mbStopped = true;
}

bool LocalMapping::isFinished()
{
    std::unique_lock<std::mutex> lock(mMutexFinish);
    return mbFinished;
}

void LocalMapping::InitializeIMU(float priorG, float priorA, bool bFIBA)
{
    if (mbResetRequested)
        return;

    float minTime;
    int nMinKF;
    if (mbMonocular)
    {
        minTime = 2.0;
        nMinKF = 10;
    }
    else
    {
        minTime = 1.0;
        nMinKF = 10;
    }


    if(mpAtlas->KeyFramesInMap()<nMinKF)
        return;

    // Retrieve all keyframe in temporal order
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while(pKF->mPrevKF)
    {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    vector<KeyFrame*> vpKF(lpKF.begin(),lpKF.end());

    if(vpKF.size()<nMinKF)
        return;

    mFirstTs=vpKF.front()->mTimeStamp;
    if(mpCurrentKeyFrame->mTimeStamp-mFirstTs<minTime)
        return;

    bInitializing = true;
    const bool map_was_imu_initialized = mpAtlas->isImuInitialized();
    const std::string alignment_event =
        !map_was_imu_initialized ? "init_accepted" : (priorA != 0.f ? "viba_1" : "viba_2");
    if (!map_was_imu_initialized)
    {
      ++mInertialInitAttemptId;
      mLastInitLifecycleCommitted = false;
    }

    while(CheckNewKeyFrames())
    {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    const int N = vpKF.size();
    IMU::Bias b(0,0,0,0,0,0);

    // Compute and KF velocities mRwg estimation
    if (!mpCurrentKeyFrame->GetMap()->isImuInitialized())
    {
        Eigen::Matrix3f Rwg;
        Eigen::Vector3f dirG;
        dirG.setZero();
        for(vector<KeyFrame*>::iterator itKF = vpKF.begin(); itKF!=vpKF.end(); itKF++)
        {
            if (!(*itKF)->mpImuPreintegrated)
                continue;
            if (!(*itKF)->mPrevKF)
                continue;

            dirG -= (*itKF)->mPrevKF->GetImuRotation() * (*itKF)->mpImuPreintegrated->GetUpdatedDeltaVelocity();
            Eigen::Vector3f _vel = ((*itKF)->GetImuPosition() - (*itKF)->mPrevKF->GetImuPosition())/(*itKF)->mpImuPreintegrated->dT;
            (*itKF)->SetVelocity(_vel);
            (*itKF)->mPrevKF->SetVelocity(_vel);
        }

        dirG = dirG/dirG.norm();
        Eigen::Vector3f gI(0.0f, 0.0f, -1.0f);
        Eigen::Vector3f v = gI.cross(dirG);
        const float nv = v.norm();
        const float cosg = gI.dot(dirG);
        const float ang = acos(cosg);
        Eigen::Vector3f vzg = v*ang/nv;
        Rwg = Sophus::SO3f::exp(vzg).matrix();
        mRwg = Rwg.cast<double>();
        mTinit = mpCurrentKeyFrame->mTimeStamp-mFirstTs;
    }
    else
    {
        mRwg = Eigen::Matrix3d::Identity();
        mbg = mpCurrentKeyFrame->GetGyroBias().cast<double>();
        mba = mpCurrentKeyFrame->GetAccBias().cast<double>();
    }

    mScale=1.0;

    mInitTime = mpTracker->mLastFrame.mTimeStamp-vpKF.front()->mTimeStamp;

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale, mbg, mba, mbMonocular, infoInertial, false, false, priorG, priorA);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (mScale<1e-1)
    {
      cout << "MI init_rejected:" << " reason=scale_too_small" << " s=" << mScale << " kfs=" << N
           << " map=" << mpAtlas->GetCurrentMap()->MapPointsInMap() << " dt=" << mInitTime
           << " ba=" << mba.norm() << " bg=" << mbg.norm() << endl;
      cout << "scale too small" << endl;
      bInitializing = false;
      return;
    }

    cout << "MI init_candidate:" << " s=" << mScale << " kfs=" << N
         << " map=" << mpAtlas->GetCurrentMap()->MapPointsInMap() << " dt=" << mInitTime
         << " ba=" << mba.norm() << " bg=" << mbg.norm() << " priorG=" << priorG
         << " priorA=" << priorA << " fiba=" << bFIBA << endl;
    if (!map_was_imu_initialized)
    {
      std::cout << "MI init_lifecycle:" << " attempt=" << mInertialInitAttemptId
                << " phase=before_align after_align=0 map_imu=0"
                << " tracking_state=" << mpTracker->mState << " bad_imu=" << (mbBadImu ? 1 : 0)
                << " motion=pending should_publish=0 committed=0 reset=0" << " mTinit=" << mTinit
                << " kfs=" << N << " map=" << mpAtlas->GetCurrentMap()->MapPointsInMap()
                << std::endl;
    }

    const Eigen::Vector3d g_old_world = mRwg * Eigen::Vector3d(0.0, 0.0, -1.0);
    const Eigen::Vector3d g_new_world(0.0, 0.0, -1.0);
    const Eigen::Vector3d align_before_gravity =
        map_was_imu_initialized ? g_new_world : g_old_world;
    const std::string align_before_gravity_frame =
        map_was_imu_initialized ? "g_new_world" : "g_old_world";
    const AlignmentTraceSnapshot align_before =
        CaptureAlignmentTraceSnapshot(mpTracker, mpCurrentKeyFrame, mScale, align_before_gravity,
                                      align_before_gravity_frame, vpKF);
    AlignmentTraceSnapshot full_ba_before;

    // Before this line we are not changing the map
    {
        std::unique_lock<std::mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
        if ((fabs(mScale - 1.f) > 0.00001) || !mbMonocular) {
            Sophus::SE3f Twg(mRwg.cast<float>().transpose(), Eigen::Vector3f::Zero());
            mpAtlas->GetCurrentMap()->ApplyScaledRotation(Twg, mScale, true);
            mpTracker->UpdateFrameIMU(mScale, vpKF[0]->GetImuBias(), mpCurrentKeyFrame);
        }

        const AlignmentTraceSnapshot align_after = CaptureAlignmentTraceSnapshot(
            mpTracker, mpCurrentKeyFrame, mScale, g_new_world, "g_new_world", vpKF);
        LogAlignmentTrace(alignment_event, align_before, align_after, mRwg, mScale);
        if (!map_was_imu_initialized)
        {
          std::cout << "MI init_lifecycle:" << " attempt=" << mInertialInitAttemptId
                    << " phase=after_align after_align=1 map_imu="
                    << (mpAtlas->isImuInitialized() ? 1 : 0)
                    << " tracking_state=" << mpTracker->mState << " bad_imu=" << (mbBadImu ? 1 : 0)
                    << " motion=pending should_publish=0 committed=0 reset=0" << " dR="
                    << RotationAngleRad(align_after.Tcw.rotationMatrix(),
                                        align_before.Tcw.rotationMatrix()) *
                           180.0 / M_PI
                    << " s=" << mScale << std::endl;
        }
        mLastAlignmentEvent = alignment_event;
        mLastAlignmentRotationDeg =
            RotationAngleRad(align_after.Tcw.rotationMatrix(), align_before.Tcw.rotationMatrix()) *
            180.0 / M_PI;
        mLastAlignmentScale = mScale;
        if (bFIBA)
          full_ba_before = align_after;

        // Check if initialization OK
        if (!mpAtlas->isImuInitialized())
            for (int i = 0; i < N; i++) {
                KeyFrame *pKF2 = vpKF[i];
                pKF2->bImu = true;
            }
    }

    mpTracker->UpdateFrameIMU(1.0,vpKF[0]->GetImuBias(),mpCurrentKeyFrame);
    if (!mpAtlas->isImuInitialized())
    {
        mpAtlas->SetImuInitialized();
        mpTracker->t0IMU = mpTracker->mCurrentFrame.mTimeStamp;
        mpCurrentKeyFrame->bImu = true;
        mInitialImuKeyframeId = mpCurrentKeyFrame->mnId;
        std::cout << "MI init_lifecycle:" << " attempt=" << mInertialInitAttemptId
                  << " phase=set_imu_initialized after_align=1 map_imu=1"
                  << " tracking_state=" << mpTracker->mState << " bad_imu=" << (mbBadImu ? 1 : 0)
                  << " motion=pending should_publish=0 committed=0 reset=0" << " mTinit=" << mTinit
                  << " keyframe=" << mpCurrentKeyFrame->mnId << std::endl;
    }

    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    if (bFIBA)
    {
        if (priorA!=0.f)
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, mpCurrentKeyFrame->mnId, NULL, true, priorG, priorA);
        else
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, mpCurrentKeyFrame->mnId, NULL, false);
    }

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    Verbose::PrintMess("Global Bundle Adjustment finished\nUpdating map ...", Verbose::VERBOSITY_NORMAL);

    // Get Map Mutex
    std::unique_lock<std::mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);

    unsigned long GBAid = mpCurrentKeyFrame->mnId;

    // Process keyframes in the queue
    while(CheckNewKeyFrames())
    {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    // Correct keyframes starting at map first keyframe
    list<KeyFrame*> lpKFtoCheck(mpAtlas->GetCurrentMap()->mvpKeyFrameOrigins.begin(),mpAtlas->GetCurrentMap()->mvpKeyFrameOrigins.end());

    while(!lpKFtoCheck.empty())
    {
        KeyFrame* pKF = lpKFtoCheck.front();
        const set<KeyFrame*> sChilds = pKF->GetChilds();
        Sophus::SE3f Twc = pKF->GetPoseInverse();
        for(set<KeyFrame*>::const_iterator sit=sChilds.begin();sit!=sChilds.end();sit++)
        {
            KeyFrame* pChild = *sit;
            if(!pChild || pChild->isBad())
                continue;

            if(pChild->mnBAGlobalForKF!=GBAid)
            {
                Sophus::SE3f Tchildc = pChild->GetPose() * Twc;
                pChild->mTcwGBA = Tchildc * pKF->mTcwGBA;

                Sophus::SO3f Rcor = pChild->mTcwGBA.so3().inverse() * pChild->GetPose().so3();
                if(pChild->isVelocitySet()){
                    pChild->mVwbGBA = Rcor * pChild->GetVelocity();
                }
                else {
                    Verbose::PrintMess("Child velocity empty!! ", Verbose::VERBOSITY_NORMAL);
                }

                pChild->mBiasGBA = pChild->GetImuBias();
                pChild->mnBAGlobalForKF = GBAid;

            }
            lpKFtoCheck.push_back(pChild);
        }

        pKF->mTcwBefGBA = pKF->GetPose();
        pKF->SetPose(pKF->mTcwGBA);

        if(pKF->bImu)
        {
            pKF->mVwbBefGBA = pKF->GetVelocity();
            pKF->SetVelocity(pKF->mVwbGBA);
            pKF->SetNewBias(pKF->mBiasGBA);
        } else {
            cout << "KF " << pKF->mnId << " not set to inertial!! \n";
        }

        lpKFtoCheck.pop_front();
    }

    // Correct MapPoints
    const vector<MapPoint*> vpMPs = mpAtlas->GetCurrentMap()->GetAllMapPoints();

    for(size_t i=0; i<vpMPs.size(); i++)
    {
        MapPoint* pMP = vpMPs[i];

        if(pMP->isBad())
            continue;

        if(pMP->mnBAGlobalForKF==GBAid)
        {
            // If optimized by Global BA, just update
            pMP->SetWorldPos(pMP->mPosGBA);
        }
        else
        {
            // Update according to the correction of its reference keyframe
            KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();

            if(pRefKF->mnBAGlobalForKF!=GBAid)
                continue;

            // Map to non-corrected camera
            Eigen::Vector3f Xc = pRefKF->mTcwBefGBA * pMP->GetWorldPos();

            // Backproject using corrected camera
            pMP->SetWorldPos(pRefKF->GetPoseInverse() * Xc);
        }
    }

    Verbose::PrintMess("Map updated!", Verbose::VERBOSITY_NORMAL);

    mnKFs=vpKF.size();
    mIdxInit++;

    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
    {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    mpTracker->mState=Tracking::OK;
    if (bFIBA)
    {
      const AlignmentTraceSnapshot full_ba_after = CaptureAlignmentTraceSnapshot(
          mpTracker, mpCurrentKeyFrame, mScale, g_new_world, "g_new_world", vpKF);
      LogAlignmentTrace("full_inertial_ba", full_ba_before, full_ba_after, mRwg, mScale);
      mLastAlignmentEvent = "full_inertial_ba";
      mLastAlignmentRotationDeg = RotationAngleRad(full_ba_after.Tcw.rotationMatrix(),
                                                   full_ba_before.Tcw.rotationMatrix()) *
                                  180.0 / M_PI;
      mLastAlignmentScale = mScale;
    }

    bInitializing = false;

    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}

void LocalMapping::ScaleRefinement()
{
    // Minimum number of keyframes to compute a solution
    // Minimum time (seconds) between first and last keyframe to compute a solution. Make the difference between monocular and stereo
    // std::unique_lock<std::mutex> lock0(mMutexImuInit);
    if (mbResetRequested)
        return;

    // Retrieve all keyframes in temporal order
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while(pKF->mPrevKF)
    {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    vector<KeyFrame*> vpKF(lpKF.begin(),lpKF.end());

    while(CheckNewKeyFrames())
    {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    const int N = vpKF.size();

    mRwg = Eigen::Matrix3d::Identity();
    mScale=1.0;

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (mScale<1e-1) // 1e-1
    {
        cout << "scale too small" << endl;
        bInitializing=false;
        return;
    }
    
    Sophus::SO3d so3wg(mRwg);
    // Before this line we are not changing the map
    std::unique_lock<std::mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    if ((fabs(mScale-1.f)>0.002)||!mbMonocular)
    {
        Sophus::SE3f Tgw(mRwg.cast<float>().transpose(),Eigen::Vector3f::Zero());
        mpAtlas->GetCurrentMap()->ApplyScaledRotation(Tgw,mScale,true);
        mpTracker->UpdateFrameIMU(mScale,mpCurrentKeyFrame->GetImuBias(),mpCurrentKeyFrame);
    }
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
    {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    double t_inertial_only = std::chrono::duration_cast<std::chrono::duration<double> >(t1 - t0).count();

    // To perform pose-inertial opt w.r.t. last keyframe
    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}



bool LocalMapping::IsInitializing()
{
    return bInitializing;
}


double LocalMapping::GetCurrKFTime()
{

    if (mpCurrentKeyFrame)
    {
        return mpCurrentKeyFrame->mTimeStamp;
    }
    else
        return 0.0;
}

KeyFrame* LocalMapping::GetCurrKF()
{
    return mpCurrentKeyFrame;
}

} //namespace ORB_SLAM
