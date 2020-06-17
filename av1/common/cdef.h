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
#ifndef AOM_AV1_COMMON_CDEF_H_
#define AOM_AV1_COMMON_CDEF_H_

#define CDEF_STRENGTH_BITS 6

#define CDEF_PRI_STRENGTHS 16
#define CDEF_SEC_STRENGTHS 4

#include "config/aom_config.h"

#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "av1/common/av1_common_int.h"
#include "av1/common/cdef_block.h"

static INLINE int sign(int i) { return i < 0 ? -1 : 1; }

static INLINE int constrain(int diff, int threshold, int damping) {
  if (!threshold) return 0;

  const int shift = AOMMAX(0, damping - get_msb(threshold));
  return sign(diff) *
         AOMMIN(abs(diff), AOMMAX(0, threshold - (abs(diff) >> shift)));
}

#ifdef __cplusplus
extern "C" {
#endif

int av1_cdef_compute_sb_list(const CommonModeInfoParams *const mi_params,
                             int mi_row, int mi_col, cdef_list *dlist,
                             BLOCK_SIZE bsize);
void av1_cdef_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm, MACROBLOCKD *xd);

/*!\brief AV1 CDEF parameter search
 *
 * \ingroup in_loop_cdef
 *
 * Searches for optimal CDEF parameters for frame
 *
 * \param[in]      frame        Compressed frame buffer
 * \param[in]      ref          Source frame buffer
 * \param[in,out]  cm           Pointer to top level common structure
 * \param[in]      xd           Pointer to common current coding block structure
 * \param[in]      pick_method  The method used to select params
 * \param[in]      rdmult       rd multiplier to use in making param choices
 *
 * \par pick_method includes:
 * \arg \c CDEF_FULL_SEARCH: Do full search
 * \arg \c CDEF_FAST_SEARCH_LVL1: Search among a subset of all possible filters.
 * \arg \c CDEF_FAST_SEARCH_LVL2: Search reduced subset of filters than Level 1.
 * \arg \c CDEF_FAST_SEARCH_LVL3: Search reduced subset of secondary filters
 *                                than Level 2.
 * \arg \c CDEF_PICK_FROM_Q: Estimate filter strength based on quantizer.
 *
 * \return Nothing is returned. Instead, optimal CDEF parameters are stored
 * in the \c cdef_info structure inside \c cm:
 * \arg \c cdef_bits: Bits of strength parameters
 * \arg \c nb_cdef_strengths: Number of strength parameters
 * \arg \c cdef_strengths: list of \c nb_cdef_strengths strength parameters
 * for the luma plane.
 * \arg \c uv_cdef_strengths: list of \c nb_cdef_strengths strength parameters
 * for the chroma planes.
 *
 */
// TODO(any): convert type of pick_method to CDEF_PICK_METHOD
void av1_cdef_search(const YV12_BUFFER_CONFIG *frame,
                     const YV12_BUFFER_CONFIG *ref, AV1_COMMON *cm,
                     MACROBLOCKD *xd, int pick_method, int rdmult);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AOM_AV1_COMMON_CDEF_H_
