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

#ifndef AOM_AOM_DSP_FLOW_ESTIMATION_DISFLOW_H_
#define AOM_AOM_DSP_FLOW_ESTIMATION_DISFLOW_H_

#include "aom_dsp/flow_estimation/flow_estimation.h"
#include "aom_scale/yv12config.h"

#ifdef __cplusplus
extern "C" {
#endif

int av1_compute_global_motion_disflow_based(
    TransformationType type, ImagePyramid *frm_pyramid, int *frm_corners,
    int num_frm_corners, ImagePyramid *ref_pyramid, int *num_inliers_by_motion,
    MotionModel *params_by_motion, int num_motions);

#ifdef __cplusplus
}
#endif

#endif  // AOM_AOM_DSP_FLOW_ESTIMATION_DISFLOW_H_
