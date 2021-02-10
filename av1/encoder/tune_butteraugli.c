/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>

#include "av1/encoder/tune_butteraugli.h"

#include "aom_dsp/butteraugli.h"
#include "aom_ports/system_state.h"
#include "av1/encoder/rdopt.h"
#include "av1/encoder/extend.h"

void av1_setup_butteraugli_recon(AV1_COMP *cpi,
                                 const YV12_BUFFER_CONFIG *recon) {
  YV12_BUFFER_CONFIG *const dst = &cpi->butteraugli_info.recon;
  AV1_COMMON *const cm = &cpi->common;
  const int width = recon->y_width;
  const int height = recon->y_height;
  if (dst->buffer_alloc_sz == 0) {
    aom_alloc_frame_buffer(
        dst, width, height, 1, 1, cm->seq_params.use_highbitdepth,
        cpi->oxcf.border_in_pixels, cm->features.byte_alignment);
  }
  av1_copy_and_extend_frame(recon, dst);
  cpi->butteraugli_info.recon_set = true;
}

void av1_set_mb_butteraugli_rdmult_scaling(AV1_COMP *cpi) {
  if (!cpi->butteraugli_info.recon_set) {
    return;
  }
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  YV12_BUFFER_CONFIG *source = cpi->source;
  const int width = source->y_width;
  const int height = source->y_height;
  const int bit_depth = cpi->td.mb.e_mbd.bd;

  aom_clear_system_state();
  const YV12_BUFFER_CONFIG *recon = &cpi->butteraugli_info.recon;
  float *diffmap;
  CHECK_MEM_ERROR(cm, diffmap, aom_malloc(width * height * sizeof(*diffmap)));
  aom_calc_butteraugli(source, recon, bit_depth, diffmap);

  const int block_size = BLOCK_8X8;
  const int num_mi_w = mi_size_wide[block_size];
  const int num_mi_h = mi_size_high[block_size];
  const int num_cols = (mi_params->mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (mi_params->mi_rows + num_mi_h - 1) / num_mi_h;
  const int block_w = num_mi_w << 2;
  const int block_h = num_mi_h << 2;
  double log_sum = 0.0;
  double blk_count = 0.0;

  // Loop through each block.
  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      const int y_start = row * block_h;
      const int x_start = col * block_w;
      float dbutteraugli = 0.0f;
      float dmse = 0.0f;

      // Loop through each pixel.
      for (int y = y_start; y < y_start + block_h && y < height; y++) {
        for (int x = x_start; x < x_start + block_w && x < width; x++) {
          dbutteraugli += powf(diffmap[y * width + x], 12.0f);
          float px_diff = source->y_buffer[y * source->y_stride + x] -
                          recon->y_buffer[y * recon->y_stride + x];
          dmse += px_diff * px_diff;
        }
      }
      for (int y = y_start; y < y_start + block_h && y < height; y += 2) {
        for (int x = x_start; x < x_start + block_w && x < width; x += 2) {
          const int px_index = y / 2 * source->uv_stride + x / 2;
          const float px_diff_u =
              source->u_buffer[px_index] - recon->u_buffer[px_index];
          const float px_diff_v =
              source->v_buffer[px_index] - recon->v_buffer[px_index];
          dmse += px_diff_u * px_diff_u + px_diff_v * px_diff_v;
        }
      }

      dbutteraugli = powf(dbutteraugli, 1.0f / 12.0f);
      dmse = dmse / (2.0f * (float)block_w * (float)block_h);
      const double K = 0.4;
      const float eps = 0.01f;
      double weight;
      if (dbutteraugli < eps || dmse < eps) {
        weight = -1.0;
      } else {
        blk_count += 1.0;
        weight = dmse / dbutteraugli + K;
        log_sum += log(weight);
      }
      cpi->butteraugli_info.rdmult_scaling_factors[index] = weight;
    }
  }
  // Geometric average of the weights.
  log_sum = exp(log_sum / blk_count);

  for (int row = 0; row < num_rows; ++row) {
    for (int col = 0; col < num_cols; ++col) {
      const int index = row * num_cols + col;
      double *weight = &cpi->butteraugli_info.rdmult_scaling_factors[index];
      if (*weight <= 0.0) {
        *weight = 1.0;
      } else {
        *weight /= log_sum;
      }
    }
  }

  aom_clear_system_state();
  aom_free(diffmap);
}

void av1_set_butteraugli_rdmult(const AV1_COMP *cpi, MACROBLOCK *x,
                                BLOCK_SIZE bsize, int mi_row, int mi_col,
                                int *rdmult) {
  assert(cpi->oxcf.tune_cfg.tuning == AOM_TUNE_BUTTERAUGLI);
  if (!cpi->butteraugli_info.recon_set) {
    return;
  }
  const AV1_COMMON *const cm = &cpi->common;

  const int bsize_base = BLOCK_8X8;
  const int num_mi_w = mi_size_wide[bsize_base];
  const int num_mi_h = mi_size_high[bsize_base];
  const int num_cols = (cm->mi_params.mi_cols + num_mi_w - 1) / num_mi_w;
  const int num_rows = (cm->mi_params.mi_rows + num_mi_h - 1) / num_mi_h;
  const int num_bcols = (mi_size_wide[bsize] + num_mi_w - 1) / num_mi_w;
  const int num_brows = (mi_size_high[bsize] + num_mi_h - 1) / num_mi_h;
  double num_of_mi = 0.0;
  double geom_mean_of_scale = 0.0;

  aom_clear_system_state();
  for (int row = mi_row / num_mi_w;
       row < num_rows && row < mi_row / num_mi_w + num_brows; ++row) {
    for (int col = mi_col / num_mi_h;
         col < num_cols && col < mi_col / num_mi_h + num_bcols; ++col) {
      const int index = row * num_cols + col;
      geom_mean_of_scale +=
          log(cpi->butteraugli_info.rdmult_scaling_factors[index]);
      num_of_mi += 1.0;
    }
  }
  geom_mean_of_scale = exp(geom_mean_of_scale / num_of_mi);

  *rdmult = (int)((double)(*rdmult) * geom_mean_of_scale + 0.5);
  *rdmult = AOMMAX(*rdmult, 0);
  av1_set_error_per_bit(&x->errorperbit, *rdmult);
  aom_clear_system_state();
}
