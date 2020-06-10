/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_ENCODER_INTRA_MODE_SEARCH_H_
#define AOM_AV1_ENCODER_INTRA_MODE_SEARCH_H_

#include "av1/encoder/encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

// Keeps track of the intra-mode search process during inter frame coding.
typedef struct IntraModeSearchState {
  // The best luma intra-mode found so far
  PREDICTION_MODE best_intra_mode;

  // Terminate intra-mode search
  int skip_intra_modes;
  // Skip the directional mode
  int dir_mode_skip_mask_ready;
  uint8_t directional_mode_skip_mask[INTRA_MODES];

  // Saving a copy of the current best chroma prediction and their statistics so
  // we don't have to recompute it every time handle_inter_mode is called.
  int rate_uv_intra;
  int rate_uv_tokenonly;
  int64_t dist_uvs;
  int skip_uvs;
  UV_PREDICTION_MODE mode_uv;
  PALETTE_MODE_INFO pmi_uv;
  int8_t uv_angle_delta;

  // Keep track of the best intra rd so it can be used in compound mode.
  int64_t best_pred_rd[REFERENCE_MODES];
} IntraModeSearchState;

// Handles an intra-mode prediction when the current frame is an inter frame.
// Its inter-mode counterpart is handle_inter_mode. This function does not
// support palette mode prediction, which uses on av1_search_palette_mode.
// The current mode being searched is defined by x->e_mbd.mi[0]->mode.
int64_t av1_handle_intra_mode(IntraModeSearchState *intra_search_state,
                              const AV1_COMP *cpi, MACROBLOCK *x,
                              BLOCK_SIZE bsize, int ref_frame_cost,
                              const PICK_MODE_CONTEXT *ctx, int disable_skip,
                              RD_STATS *rd_stats, RD_STATS *rd_stats_y,
                              RD_STATS *rd_stats_uv, int64_t best_rd,
                              int64_t *best_intra_rd);

// Handles palette-mode search when the current frame is an inter frame. In
// addition to searching palette-mode in the luma channel, this always searches
// over all possible modes for the chroma channel by calling
// av1_rd_pick_intra_sbuv_mode.
int av1_search_palette_mode(const AV1_COMP *cpi, MACROBLOCK *x,
                            RD_STATS *this_rd_cost, PICK_MODE_CONTEXT *ctx,
                            BLOCK_SIZE bsize, MB_MODE_INFO *const mbmi,
                            PALETTE_MODE_INFO *const pmi,
                            unsigned int *ref_costs_single,
                            IntraModeSearchState *intra_search_state,
                            int64_t best_rd);

// Performs intra-mode search on the luma channel when the current frame is
// intra-only. This function does not search intra-bc mode, but it does search
// palette and filter_intra.
int64_t av1_rd_pick_intra_sby_mode(const AV1_COMP *const cpi, MACROBLOCK *x,
                                   int *rate, int *rate_tokenonly,
                                   int64_t *distortion, int *skippable,
                                   BLOCK_SIZE bsize, int64_t best_rd,
                                   PICK_MODE_CONTEXT *ctx);

// Performs intra-mode search on the chroma channels. Just like its luma
// counterpart, this function searches over palette mode as well (filter_intra
// is not available for chroma channels). Unlike its luma chroma part, this
// function is used by both inter and intra frames.
int64_t av1_rd_pick_intra_sbuv_mode(const AV1_COMP *const cpi, MACROBLOCK *x,
                                    int *rate, int *rate_tokenonly,
                                    int64_t *distortion, int *skippable,
                                    BLOCK_SIZE bsize, TX_SIZE max_tx_size);

// Returns the number of colors in 'src'. This is primarily used by palette
// mode for screen content encoding.
int av1_count_colors(const uint8_t *src, int stride, int rows, int cols,
                     int *val_count);

// Same as av1_count_colors(), but for high-bitdepth mode.
int av1_count_colors_highbd(const uint8_t *src8, int stride, int rows, int cols,
                            int bit_depth, int *val_count);

// Resets palette color map for chroma channels.
void av1_restore_uv_color_map(const AV1_COMP *const cpi, MACROBLOCK *x);

#endif  // AOM_AV1_ENCODER_INTRA_MODE_SEARCH_H_
