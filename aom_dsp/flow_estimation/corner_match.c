/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "config/aom_dsp_rtcd.h"

#include "aom_dsp/flow_estimation/corner_detect.h"
#include "aom_dsp/flow_estimation/corner_match.h"
#include "aom_dsp/flow_estimation/flow_estimation.h"
#include "aom_dsp/flow_estimation/ransac.h"
#include "aom_dsp/pyramid.h"
#include "aom_scale/yv12config.h"

#define SEARCH_SZ 9
#define SEARCH_SZ_BY2 ((SEARCH_SZ - 1) / 2)

#define THRESHOLD_NCC 0.75

/* Compute var(im) * MATCH_SZ_SQ over a MATCH_SZ by MATCH_SZ window of im,
   centered at (x, y).
*/
static double compute_variance(const unsigned char *im, int stride, int x,
                               int y) {
  int sum = 0;
  int sumsq = 0;
  int var;
  int i, j;
  for (i = 0; i < MATCH_SZ; ++i)
    for (j = 0; j < MATCH_SZ; ++j) {
      sum += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
      sumsq += im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)] *
               im[(i + y - MATCH_SZ_BY2) * stride + (j + x - MATCH_SZ_BY2)];
    }
  var = sumsq * MATCH_SZ_SQ - sum * sum;
  return (double)var;
}

/* Compute corr(im1, im2) * MATCH_SZ * stddev(im1), where the
   correlation/standard deviation are taken over MATCH_SZ by MATCH_SZ windows
   of each image, centered at (x1, y1) and (x2, y2) respectively.
*/
double av1_compute_cross_correlation_c(const unsigned char *im1, int stride1,
                                       int x1, int y1, const unsigned char *im2,
                                       int stride2, int x2, int y2) {
  int v1, v2;
  int sum1 = 0;
  int sum2 = 0;
  int sumsq2 = 0;
  int cross = 0;
  int var2, cov;
  int i, j;
  for (i = 0; i < MATCH_SZ; ++i)
    for (j = 0; j < MATCH_SZ; ++j) {
      v1 = im1[(i + y1 - MATCH_SZ_BY2) * stride1 + (j + x1 - MATCH_SZ_BY2)];
      v2 = im2[(i + y2 - MATCH_SZ_BY2) * stride2 + (j + x2 - MATCH_SZ_BY2)];
      sum1 += v1;
      sum2 += v2;
      sumsq2 += v2 * v2;
      cross += v1 * v2;
    }
  var2 = sumsq2 * MATCH_SZ_SQ - sum2 * sum2;
  cov = cross * MATCH_SZ_SQ - sum1 * sum2;
  return cov / sqrt((double)var2);
}

static int is_eligible_point(int pointx, int pointy, int width, int height) {
  return (pointx >= MATCH_SZ_BY2 && pointy >= MATCH_SZ_BY2 &&
          pointx + MATCH_SZ_BY2 < width && pointy + MATCH_SZ_BY2 < height);
}

static int is_eligible_distance(int point1x, int point1y, int point2x,
                                int point2y, int width, int height) {
  const int thresh = (width < height ? height : width) >> 4;
  return ((point1x - point2x) * (point1x - point2x) +
          (point1y - point2y) * (point1y - point2y)) <= thresh * thresh;
}

static void improve_correspondence(const unsigned char *frm,
                                   const unsigned char *ref, int width,
                                   int height, int frm_stride, int ref_stride,
                                   Correspondence *correspondences,
                                   int num_correspondences) {
  int i;
  for (i = 0; i < num_correspondences; ++i) {
    int x, y, best_x = 0, best_y = 0;
    double best_match_ncc = 0.0;
    for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y) {
      for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
        double match_ncc;
        if (!is_eligible_point(correspondences[i].rx + x,
                               correspondences[i].ry + y, width, height))
          continue;
        if (!is_eligible_distance(correspondences[i].x, correspondences[i].y,
                                  correspondences[i].rx + x,
                                  correspondences[i].ry + y, width, height))
          continue;
        match_ncc = av1_compute_cross_correlation(
            frm, frm_stride, correspondences[i].x, correspondences[i].y, ref,
            ref_stride, correspondences[i].rx + x, correspondences[i].ry + y);
        if (match_ncc > best_match_ncc) {
          best_match_ncc = match_ncc;
          best_y = y;
          best_x = x;
        }
      }
    }
    correspondences[i].rx += best_x;
    correspondences[i].ry += best_y;
  }
  for (i = 0; i < num_correspondences; ++i) {
    int x, y, best_x = 0, best_y = 0;
    double best_match_ncc = 0.0;
    for (y = -SEARCH_SZ_BY2; y <= SEARCH_SZ_BY2; ++y)
      for (x = -SEARCH_SZ_BY2; x <= SEARCH_SZ_BY2; ++x) {
        double match_ncc;
        if (!is_eligible_point(correspondences[i].x + x,
                               correspondences[i].y + y, width, height))
          continue;
        if (!is_eligible_distance(
                correspondences[i].x + x, correspondences[i].y + y,
                correspondences[i].rx, correspondences[i].ry, width, height))
          continue;
        match_ncc = av1_compute_cross_correlation(
            ref, ref_stride, correspondences[i].rx, correspondences[i].ry, frm,
            frm_stride, correspondences[i].x + x, correspondences[i].y + y);
        if (match_ncc > best_match_ncc) {
          best_match_ncc = match_ncc;
          best_y = y;
          best_x = x;
        }
      }
    correspondences[i].x += best_x;
    correspondences[i].y += best_y;
  }
}

int aom_determine_correspondence(const unsigned char *src, int *src_corners,
                                 int num_src_corners, const unsigned char *ref,
                                 int *ref_corners, int num_ref_corners,
                                 int width, int height, int src_stride,
                                 int ref_stride, int *correspondence_pts) {
  // TODO(sarahparker) Improve this to include 2-way match
  int i, j;
  Correspondence *correspondences = (Correspondence *)correspondence_pts;
  int num_correspondences = 0;
  for (i = 0; i < num_src_corners; ++i) {
    double best_match_ncc = 0.0;
    double template_norm;
    int best_match_j = -1;
    if (!is_eligible_point(src_corners[2 * i], src_corners[2 * i + 1], width,
                           height))
      continue;
    for (j = 0; j < num_ref_corners; ++j) {
      double match_ncc;
      if (!is_eligible_point(ref_corners[2 * j], ref_corners[2 * j + 1], width,
                             height))
        continue;
      if (!is_eligible_distance(src_corners[2 * i], src_corners[2 * i + 1],
                                ref_corners[2 * j], ref_corners[2 * j + 1],
                                width, height))
        continue;
      match_ncc = av1_compute_cross_correlation(
          src, src_stride, src_corners[2 * i], src_corners[2 * i + 1], ref,
          ref_stride, ref_corners[2 * j], ref_corners[2 * j + 1]);
      if (match_ncc > best_match_ncc) {
        best_match_ncc = match_ncc;
        best_match_j = j;
      }
    }
    // Note: We want to test if the best correlation is >= THRESHOLD_NCC,
    // but need to account for the normalization in
    // av1_compute_cross_correlation.
    template_norm = compute_variance(src, src_stride, src_corners[2 * i],
                                     src_corners[2 * i + 1]);
    if (best_match_ncc > THRESHOLD_NCC * sqrt(template_norm)) {
      correspondences[num_correspondences].x = src_corners[2 * i];
      correspondences[num_correspondences].y = src_corners[2 * i + 1];
      correspondences[num_correspondences].rx = ref_corners[2 * best_match_j];
      correspondences[num_correspondences].ry =
          ref_corners[2 * best_match_j + 1];
      num_correspondences++;
    }
  }
  improve_correspondence(src, ref, width, height, src_stride, ref_stride,
                         correspondences, num_correspondences);
  return num_correspondences;
}

static bool get_inliers_from_indices(MotionModel *params,
                                     int *correspondences) {
  int *inliers_tmp = (int *)aom_calloc(2 * MAX_CORNERS, sizeof(*inliers_tmp));
  if (!inliers_tmp) return false;

  for (int i = 0; i < params->num_inliers; i++) {
    int index = params->inliers[i];
    inliers_tmp[2 * i] = correspondences[4 * index];
    inliers_tmp[2 * i + 1] = correspondences[4 * index + 1];
  }
  memcpy(params->inliers, inliers_tmp, sizeof(*inliers_tmp) * 2 * MAX_CORNERS);
  aom_free(inliers_tmp);
  return true;
}

int av1_compute_global_motion_feature_based(
    TransformationType type, ImagePyramid *src_pyramid, int *src_corners,
    int num_src_corners, ImagePyramid *ref_pyramid, int *num_inliers_by_motion,
    MotionModel *params_by_motion, int num_motions) {
  int i;
  int num_ref_corners;
  int num_correspondences;
  int *correspondences;
  int ref_corners[2 * MAX_CORNERS];
  RansacFunc ransac = av1_get_ransac_type(type);

  assert(aom_is_pyramid_valid(src_pyramid));
  const uint8_t *src_buffer = src_pyramid->layers[0].buffer;
  const int src_width = src_pyramid->layers[0].width;
  const int src_height = src_pyramid->layers[0].height;
  const int src_stride = src_pyramid->layers[0].stride;

  assert(aom_is_pyramid_valid(ref_pyramid));
  const uint8_t *ref_buffer = ref_pyramid->layers[0].buffer;
  const int ref_width = ref_pyramid->layers[0].width;
  const int ref_height = ref_pyramid->layers[0].height;
  const int ref_stride = ref_pyramid->layers[0].stride;

  num_ref_corners = av1_fast_corner_detect(
      ref_buffer, ref_width, ref_height, ref_stride, ref_corners, MAX_CORNERS);

  // find correspondences between the two images
  correspondences =
      (int *)aom_malloc(num_src_corners * 4 * sizeof(*correspondences));
  if (!correspondences) return 0;
  num_correspondences = aom_determine_correspondence(
      src_buffer, (int *)src_corners, num_src_corners, ref_buffer,
      (int *)ref_corners, num_ref_corners, src_width, src_height, src_stride,
      ref_stride, correspondences);

  ransac(correspondences, num_correspondences, num_inliers_by_motion,
         params_by_motion, num_motions);

  // Set num_inliers = 0 for motions with too few inliers so they are ignored.
  for (i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] < MIN_INLIER_PROB * num_correspondences ||
        num_correspondences == 0) {
      num_inliers_by_motion[i] = 0;
    } else if (!get_inliers_from_indices(&params_by_motion[i],
                                         correspondences)) {
      aom_free(correspondences);
      return 0;
    }
  }

  aom_free(correspondences);

  // Return true if any one of the motions has inliers.
  for (i = 0; i < num_motions; ++i) {
    if (num_inliers_by_motion[i] > 0) return 1;
  }
  return 0;
}
