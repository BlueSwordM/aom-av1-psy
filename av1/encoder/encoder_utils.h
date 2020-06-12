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

#ifndef AOM_AV1_ENCODER_ENCODER_UTILS_H_
#define AOM_AV1_ENCODER_ENCODER_UTILS_H_

#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define AM_SEGMENT_ID_INACTIVE 7
#define AM_SEGMENT_ID_ACTIVE 0

const int default_tx_type_probs[FRAME_UPDATE_TYPES][TX_SIZES_ALL][TX_TYPES] = {
  { { 221, 189, 214, 292, 0, 0, 0, 0, 0, 2, 38, 68, 0, 0, 0, 0 },
    { 262, 203, 216, 239, 0, 0, 0, 0, 0, 1, 37, 66, 0, 0, 0, 0 },
    { 315, 231, 239, 226, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 222, 188, 214, 287, 0, 0, 0, 0, 0, 2, 50, 61, 0, 0, 0, 0 },
    { 256, 182, 205, 282, 0, 0, 0, 0, 0, 2, 21, 76, 0, 0, 0, 0 },
    { 281, 214, 217, 222, 0, 0, 0, 0, 0, 1, 48, 41, 0, 0, 0, 0 },
    { 263, 194, 225, 225, 0, 0, 0, 0, 0, 2, 15, 100, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 170, 192, 242, 293, 0, 0, 0, 0, 0, 1, 68, 58, 0, 0, 0, 0 },
    { 199, 210, 213, 291, 0, 0, 0, 0, 0, 1, 14, 96, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 106, 69, 107, 278, 9, 15, 20, 45, 49, 23, 23, 88, 36, 74, 25, 57 },
    { 105, 72, 81, 98, 45, 49, 47, 50, 56, 72, 30, 81, 33, 95, 27, 83 },
    { 211, 105, 109, 120, 57, 62, 43, 49, 52, 58, 42, 116, 0, 0, 0, 0 },
    { 1008, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 131, 57, 98, 172, 19, 40, 37, 64, 69, 22, 41, 52, 51, 77, 35, 59 },
    { 176, 83, 93, 202, 22, 24, 28, 47, 50, 16, 12, 93, 26, 76, 17, 59 },
    { 136, 72, 89, 95, 46, 59, 47, 56, 61, 68, 35, 51, 32, 82, 26, 69 },
    { 122, 80, 87, 105, 49, 47, 46, 46, 57, 52, 13, 90, 19, 103, 15, 93 },
    { 1009, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 0, 0 },
    { 1011, 0, 0, 0, 0, 0, 0, 0, 0, 13, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 202, 20, 84, 114, 14, 60, 41, 79, 99, 21, 41, 15, 50, 84, 34, 66 },
    { 196, 44, 23, 72, 30, 22, 28, 57, 67, 13, 4, 165, 15, 148, 9, 131 },
    { 882, 0, 0, 0, 0, 0, 0, 0, 0, 142, 0, 0, 0, 0, 0, 0 },
    { 840, 0, 0, 0, 0, 0, 0, 0, 0, 184, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 } },
  { { 213, 110, 141, 269, 12, 16, 15, 19, 21, 11, 38, 68, 22, 29, 16, 24 },
    { 216, 119, 128, 143, 38, 41, 26, 30, 31, 30, 42, 70, 23, 36, 19, 32 },
    { 367, 149, 154, 154, 38, 35, 17, 21, 21, 10, 22, 36, 0, 0, 0, 0 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 219, 96, 127, 191, 21, 40, 25, 32, 34, 18, 45, 45, 33, 39, 26, 33 },
    { 296, 99, 122, 198, 23, 21, 19, 24, 25, 13, 20, 64, 23, 32, 18, 27 },
    { 275, 128, 142, 143, 35, 48, 23, 30, 29, 18, 42, 36, 18, 23, 14, 20 },
    { 239, 132, 166, 175, 36, 27, 19, 21, 24, 14, 13, 85, 9, 31, 8, 25 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1022, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 309, 25, 79, 59, 25, 80, 34, 53, 61, 25, 49, 23, 43, 64, 36, 59 },
    { 270, 57, 40, 54, 50, 42, 41, 53, 56, 28, 17, 81, 45, 86, 34, 70 },
    { 1005, 0, 0, 0, 0, 0, 0, 0, 0, 19, 0, 0, 0, 0, 0, 0 },
    { 992, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 133, 63, 55, 83, 57, 87, 58, 72, 68, 16, 24, 35, 29, 105, 25, 114 },
    { 131, 75, 74, 60, 71, 77, 65, 66, 73, 33, 21, 79, 20, 83, 18, 78 },
    { 276, 95, 82, 58, 86, 93, 63, 60, 64, 17, 38, 92, 0, 0, 0, 0 },
    { 1006, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 147, 49, 75, 78, 50, 97, 60, 67, 76, 17, 42, 35, 31, 93, 27, 80 },
    { 157, 49, 58, 75, 61, 52, 56, 67, 69, 12, 15, 79, 24, 119, 11, 120 },
    { 178, 69, 83, 77, 69, 85, 72, 77, 77, 20, 35, 40, 25, 48, 23, 46 },
    { 174, 55, 64, 57, 73, 68, 62, 61, 75, 15, 12, 90, 17, 99, 16, 86 },
    { 1008, 0, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0 },
    { 1018, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 266, 31, 63, 64, 21, 52, 39, 54, 63, 30, 52, 31, 48, 89, 46, 75 },
    { 272, 26, 32, 44, 29, 31, 32, 53, 51, 13, 13, 88, 22, 153, 16, 149 },
    { 923, 0, 0, 0, 0, 0, 0, 0, 0, 101, 0, 0, 0, 0, 0, 0 },
    { 969, 0, 0, 0, 0, 0, 0, 0, 0, 55, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
  { { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 },
    { 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 } },
  { { 158, 92, 125, 298, 12, 15, 20, 29, 31, 12, 29, 67, 34, 44, 23, 35 },
    { 147, 94, 103, 123, 45, 48, 38, 41, 46, 48, 37, 78, 33, 63, 27, 53 },
    { 268, 126, 125, 136, 54, 53, 31, 38, 38, 33, 35, 87, 0, 0, 0, 0 },
    { 1018, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 159, 72, 103, 194, 20, 35, 37, 50, 56, 21, 39, 40, 51, 61, 38, 48 },
    { 259, 86, 95, 188, 32, 20, 25, 34, 37, 13, 12, 85, 25, 53, 17, 43 },
    { 189, 99, 113, 123, 45, 59, 37, 46, 48, 44, 39, 41, 31, 47, 26, 37 },
    { 175, 110, 113, 128, 58, 38, 33, 33, 43, 29, 13, 100, 14, 68, 12, 57 },
    { 1017, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0 },
    { 1019, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 208, 22, 84, 101, 21, 59, 44, 70, 90, 25, 59, 13, 64, 67, 49, 48 },
    { 277, 52, 32, 63, 43, 26, 33, 48, 54, 11, 6, 130, 18, 119, 11, 101 },
    { 963, 0, 0, 0, 0, 0, 0, 0, 0, 61, 0, 0, 0, 0, 0, 0 },
    { 979, 0, 0, 0, 0, 0, 0, 0, 0, 45, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } }
};

const int default_obmc_probs[FRAME_UPDATE_TYPES][BLOCK_SIZES_ALL] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  106, 90, 90, 97, 67, 59, 70, 28,
    30, 38, 16, 16,  16, 0,  0,  44, 50, 26, 25 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  98, 93, 97, 68, 82, 85, 33, 30,
    33, 16, 16, 16, 16, 0,  0,  43, 37, 26, 16 },
  { 0,  0,  0,  91, 80, 76, 78, 55, 49, 24, 16,
    16, 16, 16, 16, 16, 0,  0,  29, 45, 16, 38 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0,  0,  0,  103, 89, 89, 89, 62, 63, 76, 34,
    35, 32, 19, 16,  16, 0,  0,  49, 55, 29, 19 }
};

const int default_warped_probs[FRAME_UPDATE_TYPES] = { 64, 64, 64, 64,
                                                       64, 64, 64 };

// TODO(yunqing): the default probs can be trained later from better
// performance.
const int default_switchable_interp_probs[FRAME_UPDATE_TYPES]
                                         [SWITCHABLE_FILTER_CONTEXTS]
                                         [SWITCHABLE_FILTERS] = {
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } },
                                           { { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 },
                                             { 512, 512, 512 } }
                                         };

// Mark all inactive blocks as active. Other segmentation features may be set
// so memset cannot be used, instead only inactive blocks should be reset.
static AOM_INLINE void suppress_active_map(AV1_COMP *cpi) {
  unsigned char *const seg_map = cpi->enc_seg.map;
  int i;
  if (cpi->active_map.enabled || cpi->active_map.update)
    for (i = 0;
         i < cpi->common.mi_params.mi_rows * cpi->common.mi_params.mi_cols; ++i)
      if (seg_map[i] == AM_SEGMENT_ID_INACTIVE)
        seg_map[i] = AM_SEGMENT_ID_ACTIVE;
}

static AOM_INLINE void set_mb_mi(CommonModeInfoParams *mi_params, int width,
                                 int height) {
  // Ensure that the decoded width and height are both multiples of
  // 8 luma pixels (note: this may only be a multiple of 4 chroma pixels if
  // subsampling is used).
  // This simplifies the implementation of various experiments,
  // eg. cdef, which operates on units of 8x8 luma pixels.
  const int aligned_width = ALIGN_POWER_OF_TWO(width, 3);
  const int aligned_height = ALIGN_POWER_OF_TWO(height, 3);

  mi_params->mi_cols = aligned_width >> MI_SIZE_LOG2;
  mi_params->mi_rows = aligned_height >> MI_SIZE_LOG2;
  mi_params->mi_stride = calc_mi_size(mi_params->mi_cols);

  mi_params->mb_cols = (mi_params->mi_cols + 2) >> 2;
  mi_params->mb_rows = (mi_params->mi_rows + 2) >> 2;
  mi_params->MBs = mi_params->mb_rows * mi_params->mb_cols;

  const int mi_alloc_size_1d = mi_size_wide[mi_params->mi_alloc_bsize];
  mi_params->mi_alloc_stride =
      (mi_params->mi_stride + mi_alloc_size_1d - 1) / mi_alloc_size_1d;

  assert(mi_size_wide[mi_params->mi_alloc_bsize] ==
         mi_size_high[mi_params->mi_alloc_bsize]);

#if CONFIG_LPF_MASK
  av1_alloc_loop_filter_mask(mi_params);
#endif
}

static AOM_INLINE void enc_free_mi(CommonModeInfoParams *mi_params) {
  aom_free(mi_params->mi_alloc);
  mi_params->mi_alloc = NULL;
  aom_free(mi_params->mi_grid_base);
  mi_params->mi_grid_base = NULL;
  mi_params->mi_alloc_size = 0;
  aom_free(mi_params->tx_type_map);
  mi_params->tx_type_map = NULL;
}

static AOM_INLINE void enc_set_mb_mi(CommonModeInfoParams *mi_params, int width,
                                     int height) {
  const int is_4k_or_larger = AOMMIN(width, height) >= 2160;
  mi_params->mi_alloc_bsize = is_4k_or_larger ? BLOCK_8X8 : BLOCK_4X4;

  set_mb_mi(mi_params, width, height);
}

static AOM_INLINE void stat_stage_set_mb_mi(CommonModeInfoParams *mi_params,
                                            int width, int height) {
  mi_params->mi_alloc_bsize = BLOCK_16X16;

  set_mb_mi(mi_params, width, height);
}

static AOM_INLINE void enc_setup_mi(CommonModeInfoParams *mi_params) {
  const int mi_grid_size =
      mi_params->mi_stride * calc_mi_size(mi_params->mi_rows);
  memset(mi_params->mi_alloc, 0,
         mi_params->mi_alloc_size * sizeof(*mi_params->mi_alloc));
  memset(mi_params->mi_grid_base, 0,
         mi_grid_size * sizeof(*mi_params->mi_grid_base));
  memset(mi_params->tx_type_map, 0,
         mi_grid_size * sizeof(*mi_params->tx_type_map));
}

static AOM_INLINE void dealloc_context_buffers_ext(
    MBMIExtFrameBufferInfo *mbmi_ext_info) {
  if (mbmi_ext_info->frame_base) {
    aom_free(mbmi_ext_info->frame_base);
    mbmi_ext_info->frame_base = NULL;
    mbmi_ext_info->alloc_size = 0;
  }
}

static AOM_INLINE void alloc_context_buffers_ext(
    AV1_COMMON *cm, MBMIExtFrameBufferInfo *mbmi_ext_info) {
  const CommonModeInfoParams *const mi_params = &cm->mi_params;

  const int mi_alloc_size_1d = mi_size_wide[mi_params->mi_alloc_bsize];
  const int mi_alloc_rows =
      (mi_params->mi_rows + mi_alloc_size_1d - 1) / mi_alloc_size_1d;
  const int mi_alloc_cols =
      (mi_params->mi_cols + mi_alloc_size_1d - 1) / mi_alloc_size_1d;
  const int new_ext_mi_size = mi_alloc_rows * mi_alloc_cols;

  if (new_ext_mi_size > mbmi_ext_info->alloc_size) {
    dealloc_context_buffers_ext(mbmi_ext_info);
    CHECK_MEM_ERROR(
        cm, mbmi_ext_info->frame_base,
        aom_calloc(new_ext_mi_size, sizeof(*mbmi_ext_info->frame_base)));
    mbmi_ext_info->alloc_size = new_ext_mi_size;
  }
  // The stride needs to be updated regardless of whether new allocation
  // happened or not.
  mbmi_ext_info->stride = mi_alloc_cols;
}

static AOM_INLINE void alloc_altref_frame_buffer(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;
  const SequenceHeader *const seq_params = &cm->seq_params;
  const AV1EncoderConfig *oxcf = &cpi->oxcf;

  // TODO(agrange) Check if ARF is enabled and skip allocation if not.
  if (aom_realloc_frame_buffer(
          &cpi->alt_ref_buffer, oxcf->width, oxcf->height,
          seq_params->subsampling_x, seq_params->subsampling_y,
          seq_params->use_highbitdepth, cpi->oxcf.border_in_pixels,
          cm->features.byte_alignment, NULL, NULL, NULL))
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate altref buffer");
}

static AOM_INLINE void alloc_util_frame_buffers(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const SequenceHeader *const seq_params = &cm->seq_params;
  const int byte_alignment = cm->features.byte_alignment;
  if (aom_realloc_frame_buffer(
          &cpi->last_frame_uf, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, byte_alignment, NULL, NULL, NULL))
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate last frame buffer");

  if (aom_realloc_frame_buffer(
          &cpi->trial_frame_rst, cm->superres_upscaled_width,
          cm->superres_upscaled_height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          AOM_RESTORATION_FRAME_BORDER, byte_alignment, NULL, NULL, NULL))
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate trial restored frame buffer");

  if (aom_realloc_frame_buffer(
          &cpi->scaled_source, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, byte_alignment, NULL, NULL, NULL))
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate scaled source buffer");

  if (aom_realloc_frame_buffer(
          &cpi->scaled_last_source, cm->width, cm->height,
          seq_params->subsampling_x, seq_params->subsampling_y,
          seq_params->use_highbitdepth, cpi->oxcf.border_in_pixels,
          byte_alignment, NULL, NULL, NULL))
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate scaled last source buffer");
}

static AOM_INLINE void alloc_compressor_data(AV1_COMP *cpi) {
  AV1_COMMON *cm = &cpi->common;
  TokenInfo *token_info = &cpi->token_info;

  if (av1_alloc_context_buffers(cm, cm->width, cm->height)) {
    aom_internal_error(&cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate context buffers");
  }

  if (!is_stat_generation_stage(cpi)) {
    av1_alloc_txb_buf(cpi);

    alloc_context_buffers_ext(cm, &cpi->mbmi_ext_info);
  }

  free_token_info(token_info);

  if (!is_stat_generation_stage(cpi)) {
    alloc_token_info(cm, token_info);
  }

  av1_setup_shared_coeff_buffer(&cpi->common, &cpi->td.shared_coeff_buf);
  av1_setup_sms_tree(cpi, &cpi->td);
  cpi->td.firstpass_ctx =
      av1_alloc_pmc(cm, BLOCK_16X16, &cpi->td.shared_coeff_buf);
}

static AOM_INLINE void dealloc_compressor_data(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  TokenInfo *token_info = &cpi->token_info;

  dealloc_context_buffers_ext(&cpi->mbmi_ext_info);

  aom_free(cpi->tile_data);
  cpi->tile_data = NULL;

  // Delete sementation map
  aom_free(cpi->enc_seg.map);
  cpi->enc_seg.map = NULL;

  av1_cyclic_refresh_free(cpi->cyclic_refresh);
  cpi->cyclic_refresh = NULL;

  aom_free(cpi->active_map.map);
  cpi->active_map.map = NULL;

  aom_free(cpi->ssim_rdmult_scaling_factors);
  cpi->ssim_rdmult_scaling_factors = NULL;

  aom_free(cpi->tpl_rdmult_scaling_factors);
  cpi->tpl_rdmult_scaling_factors = NULL;

  aom_free(cpi->tpl_sb_rdmult_scaling_factors);
  cpi->tpl_sb_rdmult_scaling_factors = NULL;

#if CONFIG_TUNE_VMAF
  aom_free(cpi->vmaf_info.rdmult_scaling_factors);
  cpi->vmaf_info.rdmult_scaling_factors = NULL;
#endif

  av1_release_obmc_buffers(&cpi->td.mb.obmc_buffer);

  aom_free(cpi->td.mb.inter_modes_info);
  cpi->td.mb.inter_modes_info = NULL;

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
      aom_free(cpi->td.mb.intrabc_hash_info.hash_value_buffer[i][j]);
      cpi->td.mb.intrabc_hash_info.hash_value_buffer[i][j] = NULL;
    }

  aom_free(cm->tpl_mvs);
  cm->tpl_mvs = NULL;

  aom_free(cpi->td.mb.mbmi_ext);
  cpi->td.mb.mbmi_ext = NULL;

  if (cpi->td.vt64x64) {
    aom_free(cpi->td.vt64x64);
    cpi->td.vt64x64 = NULL;
  }

  av1_free_pmc(cpi->td.firstpass_ctx, av1_num_planes(cm));
  cpi->td.firstpass_ctx = NULL;

  av1_free_ref_frame_buffers(cm->buffer_pool);
  av1_free_txb_buf(cpi);
  av1_free_context_buffers(cm);

  aom_free_frame_buffer(&cpi->last_frame_uf);
  av1_free_restoration_buffers(cm);
  aom_free_frame_buffer(&cpi->trial_frame_rst);
  aom_free_frame_buffer(&cpi->scaled_source);
  aom_free_frame_buffer(&cpi->scaled_last_source);
  aom_free_frame_buffer(&cpi->alt_ref_buffer);
  av1_lookahead_destroy(cpi->lookahead);

  free_token_info(token_info);

  av1_free_shared_coeff_buffer(&cpi->td.shared_coeff_buf);
  av1_free_sms_tree(&cpi->td);

  aom_free(cpi->td.mb.palette_buffer);
  av1_release_compound_type_rd_buffers(&cpi->td.mb.comp_rd_buffer);
  aom_free(cpi->td.mb.tmp_conv_dst);
  for (int j = 0; j < 2; ++j) {
    aom_free(cpi->td.mb.tmp_pred_bufs[j]);
  }

#if CONFIG_DENOISE
  if (cpi->denoise_and_model) {
    aom_denoise_and_model_free(cpi->denoise_and_model);
    cpi->denoise_and_model = NULL;
  }
#endif
  if (cpi->film_grain_table) {
    aom_film_grain_table_free(cpi->film_grain_table);
    cpi->film_grain_table = NULL;
  }

  for (int i = 0; i < MAX_NUM_OPERATING_POINTS; ++i) {
    aom_free(cpi->level_params.level_info[i]);
  }

  if (cpi->use_svc) av1_free_svc_cyclic_refresh(cpi);
}

static AOM_INLINE void init_buffer_indices(
    ForceIntegerMVInfo *const force_intpel_info, int *const remapped_ref_idx) {
  int fb_idx;
  for (fb_idx = 0; fb_idx < REF_FRAMES; ++fb_idx)
    remapped_ref_idx[fb_idx] = fb_idx;
  force_intpel_info->rate_index = 0;
  force_intpel_info->rate_size = 0;
}

#define HIGHBD_BFP(BT, SDF, SDAF, VF, SVF, SVAF, SDX4DF, JSDAF, JSVAF) \
  cpi->fn_ptr[BT].sdf = SDF;                                           \
  cpi->fn_ptr[BT].sdaf = SDAF;                                         \
  cpi->fn_ptr[BT].vf = VF;                                             \
  cpi->fn_ptr[BT].svf = SVF;                                           \
  cpi->fn_ptr[BT].svaf = SVAF;                                         \
  cpi->fn_ptr[BT].sdx4df = SDX4DF;                                     \
  cpi->fn_ptr[BT].jsdaf = JSDAF;                                       \
  cpi->fn_ptr[BT].jsvaf = JSVAF;

#define HIGHBD_BFP_WRAPPER(WIDTH, HEIGHT, BD)                                \
  HIGHBD_BFP(                                                                \
      BLOCK_##WIDTH##X##HEIGHT, aom_highbd_sad##WIDTH##x##HEIGHT##_bits##BD, \
      aom_highbd_sad##WIDTH##x##HEIGHT##_avg_bits##BD,                       \
      aom_highbd_##BD##_variance##WIDTH##x##HEIGHT,                          \
      aom_highbd_##BD##_sub_pixel_variance##WIDTH##x##HEIGHT,                \
      aom_highbd_##BD##_sub_pixel_avg_variance##WIDTH##x##HEIGHT,            \
      aom_highbd_sad##WIDTH##x##HEIGHT##x4d_bits##BD,                        \
      aom_highbd_dist_wtd_sad##WIDTH##x##HEIGHT##_avg_bits##BD,              \
      aom_highbd_##BD##_dist_wtd_sub_pixel_avg_variance##WIDTH##x##HEIGHT)

#define MAKE_BFP_SAD_WRAPPER(fnname)                                           \
  static unsigned int fnname##_bits8(const uint8_t *src_ptr,                   \
                                     int source_stride,                        \
                                     const uint8_t *ref_ptr, int ref_stride) { \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride);                \
  }                                                                            \
  static unsigned int fnname##_bits10(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride) {                                                        \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride) >> 2;           \
  }                                                                            \
  static unsigned int fnname##_bits12(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride) {                                                        \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride) >> 4;           \
  }

#define MAKE_BFP_SADAVG_WRAPPER(fnname)                                        \
  static unsigned int fnname##_bits8(                                          \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred);   \
  }                                                                            \
  static unsigned int fnname##_bits10(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred) >> \
           2;                                                                  \
  }                                                                            \
  static unsigned int fnname##_bits12(                                         \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,       \
      int ref_stride, const uint8_t *second_pred) {                            \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred) >> \
           4;                                                                  \
  }

#define MAKE_BFP_SAD4D_WRAPPER(fnname)                                        \
  static void fnname##_bits8(const uint8_t *src_ptr, int source_stride,       \
                             const uint8_t *const ref_ptr[], int ref_stride,  \
                             unsigned int *sad_array) {                       \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
  }                                                                           \
  static void fnname##_bits10(const uint8_t *src_ptr, int source_stride,      \
                              const uint8_t *const ref_ptr[], int ref_stride, \
                              unsigned int *sad_array) {                      \
    int i;                                                                    \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
    for (i = 0; i < 4; i++) sad_array[i] >>= 2;                               \
  }                                                                           \
  static void fnname##_bits12(const uint8_t *src_ptr, int source_stride,      \
                              const uint8_t *const ref_ptr[], int ref_stride, \
                              unsigned int *sad_array) {                      \
    int i;                                                                    \
    fnname(src_ptr, source_stride, ref_ptr, ref_stride, sad_array);           \
    for (i = 0; i < 4; i++) sad_array[i] >>= 4;                               \
  }

#define MAKE_BFP_JSADAVG_WRAPPER(fnname)                                    \
  static unsigned int fnname##_bits8(                                       \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param);                                               \
  }                                                                         \
  static unsigned int fnname##_bits10(                                      \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param) >>                                             \
           2;                                                               \
  }                                                                         \
  static unsigned int fnname##_bits12(                                      \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr,    \
      int ref_stride, const uint8_t *second_pred,                           \
      const DIST_WTD_COMP_PARAMS *jcp_param) {                              \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride, second_pred, \
                  jcp_param) >>                                             \
           4;                                                               \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad128x128)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad128x128_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad128x128x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad128x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad128x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad128x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x128)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x128_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x128x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x4x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x4x4d)

MAKE_BFP_SAD_WRAPPER(aom_highbd_sad4x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad4x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad4x16x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x4)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x4_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x4x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad8x32)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad8x32_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad8x32x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad32x8)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad32x8_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad32x8x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad16x64)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad16x64_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad16x64x4d)
MAKE_BFP_SAD_WRAPPER(aom_highbd_sad64x16)
MAKE_BFP_SADAVG_WRAPPER(aom_highbd_sad64x16_avg)
MAKE_BFP_SAD4D_WRAPPER(aom_highbd_sad64x16x4d)

MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad128x128_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad128x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x128_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad4x16_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x4_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad8x32_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad32x8_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad16x64_avg)
MAKE_BFP_JSADAVG_WRAPPER(aom_highbd_dist_wtd_sad64x16_avg)
#endif  // CONFIG_AV1_HIGHBITDEPTH

#define HIGHBD_MBFP(BT, MCSDF, MCSVF) \
  cpi->fn_ptr[BT].msdf = MCSDF;       \
  cpi->fn_ptr[BT].msvf = MCSVF;

#define HIGHBD_MBFP_WRAPPER(WIDTH, HEIGHT, BD)                    \
  HIGHBD_MBFP(BLOCK_##WIDTH##X##HEIGHT,                           \
              aom_highbd_masked_sad##WIDTH##x##HEIGHT##_bits##BD, \
              aom_highbd_##BD##_masked_sub_pixel_variance##WIDTH##x##HEIGHT)

#define MAKE_MBFP_COMPOUND_SAD_WRAPPER(fnname)                           \
  static unsigned int fnname##_bits8(                                    \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask);            \
  }                                                                      \
  static unsigned int fnname##_bits10(                                   \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask) >>          \
           2;                                                            \
  }                                                                      \
  static unsigned int fnname##_bits12(                                   \
      const uint8_t *src_ptr, int source_stride, const uint8_t *ref_ptr, \
      int ref_stride, const uint8_t *second_pred_ptr, const uint8_t *m,  \
      int m_stride, int invert_mask) {                                   \
    return fnname(src_ptr, source_stride, ref_ptr, ref_stride,           \
                  second_pred_ptr, m, m_stride, invert_mask) >>          \
           4;                                                            \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad128x128)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad128x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x128)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad4x16)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x4)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad8x32)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad32x8)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad16x64)
MAKE_MBFP_COMPOUND_SAD_WRAPPER(aom_highbd_masked_sad64x16)
#endif

#define HIGHBD_OBFP(BT, OSDF, OVF, OSVF) \
  cpi->fn_ptr[BT].osdf = OSDF;           \
  cpi->fn_ptr[BT].ovf = OVF;             \
  cpi->fn_ptr[BT].osvf = OSVF;

#define HIGHBD_OBFP_WRAPPER(WIDTH, HEIGHT, BD)                   \
  HIGHBD_OBFP(BLOCK_##WIDTH##X##HEIGHT,                          \
              aom_highbd_obmc_sad##WIDTH##x##HEIGHT##_bits##BD,  \
              aom_highbd_##BD##_obmc_variance##WIDTH##x##HEIGHT, \
              aom_highbd_##BD##_obmc_sub_pixel_variance##WIDTH##x##HEIGHT)

#define LOWBD_OBFP_WRAPPER(WIDTH, HEIGHT)                    \
  HIGHBD_OBFP(BLOCK_##WIDTH##X##HEIGHT,                      \
              aom_highbd_obmc_sad##WIDTH##x##HEIGHT##_bits8, \
              aom_highbd_obmc_variance##WIDTH##x##HEIGHT,    \
              aom_highbd_obmc_sub_pixel_variance##WIDTH##x##HEIGHT)

#define MAKE_OBFP_SAD_WRAPPER(fnname)                                     \
  static unsigned int fnname##_bits8(const uint8_t *ref, int ref_stride,  \
                                     const int32_t *wsrc,                 \
                                     const int32_t *msk) {                \
    return fnname(ref, ref_stride, wsrc, msk);                            \
  }                                                                       \
  static unsigned int fnname##_bits10(const uint8_t *ref, int ref_stride, \
                                      const int32_t *wsrc,                \
                                      const int32_t *msk) {               \
    return fnname(ref, ref_stride, wsrc, msk) >> 2;                       \
  }                                                                       \
  static unsigned int fnname##_bits12(const uint8_t *ref, int ref_stride, \
                                      const int32_t *wsrc,                \
                                      const int32_t *msk) {               \
    return fnname(ref, ref_stride, wsrc, msk) >> 4;                       \
  }

#if CONFIG_AV1_HIGHBITDEPTH
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad128x128)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad128x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x128)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad4x16)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x4)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad8x32)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad32x8)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad16x64)
MAKE_OBFP_SAD_WRAPPER(aom_highbd_obmc_sad64x16)

static AOM_INLINE void highbd_set_var_fns(AV1_COMP *const cpi) {
  AV1_COMMON *const cm = &cpi->common;
  if (cm->seq_params.use_highbitdepth) {
    switch (cm->seq_params.bit_depth) {
      case AOM_BITS_8:
        HIGHBD_BFP_WRAPPER(64, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 64, 8)
        HIGHBD_BFP_WRAPPER(32, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 32, 8)
        HIGHBD_BFP_WRAPPER(16, 4, 8)
        HIGHBD_BFP_WRAPPER(4, 16, 8)
        HIGHBD_BFP_WRAPPER(32, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 32, 8)
        HIGHBD_BFP_WRAPPER(64, 32, 8)
        HIGHBD_BFP_WRAPPER(32, 64, 8)
        HIGHBD_BFP_WRAPPER(32, 32, 8)
        HIGHBD_BFP_WRAPPER(64, 64, 8)
        HIGHBD_BFP_WRAPPER(16, 16, 8)
        HIGHBD_BFP_WRAPPER(16, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 16, 8)
        HIGHBD_BFP_WRAPPER(8, 8, 8)
        HIGHBD_BFP_WRAPPER(8, 4, 8)
        HIGHBD_BFP_WRAPPER(4, 8, 8)
        HIGHBD_BFP_WRAPPER(4, 4, 8)
        HIGHBD_BFP_WRAPPER(128, 128, 8)
        HIGHBD_BFP_WRAPPER(128, 64, 8)
        HIGHBD_BFP_WRAPPER(64, 128, 8)

        HIGHBD_MBFP_WRAPPER(128, 128, 8)
        HIGHBD_MBFP_WRAPPER(128, 64, 8)
        HIGHBD_MBFP_WRAPPER(64, 128, 8)
        HIGHBD_MBFP_WRAPPER(64, 64, 8)
        HIGHBD_MBFP_WRAPPER(64, 32, 8)
        HIGHBD_MBFP_WRAPPER(32, 64, 8)
        HIGHBD_MBFP_WRAPPER(32, 32, 8)
        HIGHBD_MBFP_WRAPPER(32, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 32, 8)
        HIGHBD_MBFP_WRAPPER(16, 16, 8)
        HIGHBD_MBFP_WRAPPER(8, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 8, 8)
        HIGHBD_MBFP_WRAPPER(4, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 4, 8)
        HIGHBD_MBFP_WRAPPER(4, 4, 8)
        HIGHBD_MBFP_WRAPPER(64, 16, 8)
        HIGHBD_MBFP_WRAPPER(16, 64, 8)
        HIGHBD_MBFP_WRAPPER(32, 8, 8)
        HIGHBD_MBFP_WRAPPER(8, 32, 8)
        HIGHBD_MBFP_WRAPPER(16, 4, 8)
        HIGHBD_MBFP_WRAPPER(4, 16, 8)

        LOWBD_OBFP_WRAPPER(128, 128)
        LOWBD_OBFP_WRAPPER(128, 64)
        LOWBD_OBFP_WRAPPER(64, 128)
        LOWBD_OBFP_WRAPPER(64, 64)
        LOWBD_OBFP_WRAPPER(64, 32)
        LOWBD_OBFP_WRAPPER(32, 64)
        LOWBD_OBFP_WRAPPER(32, 32)
        LOWBD_OBFP_WRAPPER(32, 16)
        LOWBD_OBFP_WRAPPER(16, 32)
        LOWBD_OBFP_WRAPPER(16, 16)
        LOWBD_OBFP_WRAPPER(8, 16)
        LOWBD_OBFP_WRAPPER(16, 8)
        LOWBD_OBFP_WRAPPER(8, 8)
        LOWBD_OBFP_WRAPPER(4, 8)
        LOWBD_OBFP_WRAPPER(8, 4)
        LOWBD_OBFP_WRAPPER(4, 4)
        LOWBD_OBFP_WRAPPER(64, 16)
        LOWBD_OBFP_WRAPPER(16, 64)
        LOWBD_OBFP_WRAPPER(32, 8)
        LOWBD_OBFP_WRAPPER(8, 32)
        LOWBD_OBFP_WRAPPER(16, 4)
        LOWBD_OBFP_WRAPPER(4, 16)
        break;

      case AOM_BITS_10:
        HIGHBD_BFP_WRAPPER(64, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 64, 10)
        HIGHBD_BFP_WRAPPER(32, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 32, 10)
        HIGHBD_BFP_WRAPPER(16, 4, 10)
        HIGHBD_BFP_WRAPPER(4, 16, 10)
        HIGHBD_BFP_WRAPPER(32, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 32, 10)
        HIGHBD_BFP_WRAPPER(64, 32, 10)
        HIGHBD_BFP_WRAPPER(32, 64, 10)
        HIGHBD_BFP_WRAPPER(32, 32, 10)
        HIGHBD_BFP_WRAPPER(64, 64, 10)
        HIGHBD_BFP_WRAPPER(16, 16, 10)
        HIGHBD_BFP_WRAPPER(16, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 16, 10)
        HIGHBD_BFP_WRAPPER(8, 8, 10)
        HIGHBD_BFP_WRAPPER(8, 4, 10)
        HIGHBD_BFP_WRAPPER(4, 8, 10)
        HIGHBD_BFP_WRAPPER(4, 4, 10)
        HIGHBD_BFP_WRAPPER(128, 128, 10)
        HIGHBD_BFP_WRAPPER(128, 64, 10)
        HIGHBD_BFP_WRAPPER(64, 128, 10)

        HIGHBD_MBFP_WRAPPER(128, 128, 10)
        HIGHBD_MBFP_WRAPPER(128, 64, 10)
        HIGHBD_MBFP_WRAPPER(64, 128, 10)
        HIGHBD_MBFP_WRAPPER(64, 64, 10)
        HIGHBD_MBFP_WRAPPER(64, 32, 10)
        HIGHBD_MBFP_WRAPPER(32, 64, 10)
        HIGHBD_MBFP_WRAPPER(32, 32, 10)
        HIGHBD_MBFP_WRAPPER(32, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 32, 10)
        HIGHBD_MBFP_WRAPPER(16, 16, 10)
        HIGHBD_MBFP_WRAPPER(8, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 8, 10)
        HIGHBD_MBFP_WRAPPER(4, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 4, 10)
        HIGHBD_MBFP_WRAPPER(4, 4, 10)
        HIGHBD_MBFP_WRAPPER(64, 16, 10)
        HIGHBD_MBFP_WRAPPER(16, 64, 10)
        HIGHBD_MBFP_WRAPPER(32, 8, 10)
        HIGHBD_MBFP_WRAPPER(8, 32, 10)
        HIGHBD_MBFP_WRAPPER(16, 4, 10)
        HIGHBD_MBFP_WRAPPER(4, 16, 10)

        HIGHBD_OBFP_WRAPPER(128, 128, 10)
        HIGHBD_OBFP_WRAPPER(128, 64, 10)
        HIGHBD_OBFP_WRAPPER(64, 128, 10)
        HIGHBD_OBFP_WRAPPER(64, 64, 10)
        HIGHBD_OBFP_WRAPPER(64, 32, 10)
        HIGHBD_OBFP_WRAPPER(32, 64, 10)
        HIGHBD_OBFP_WRAPPER(32, 32, 10)
        HIGHBD_OBFP_WRAPPER(32, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 32, 10)
        HIGHBD_OBFP_WRAPPER(16, 16, 10)
        HIGHBD_OBFP_WRAPPER(8, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 8, 10)
        HIGHBD_OBFP_WRAPPER(4, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 4, 10)
        HIGHBD_OBFP_WRAPPER(4, 4, 10)
        HIGHBD_OBFP_WRAPPER(64, 16, 10)
        HIGHBD_OBFP_WRAPPER(16, 64, 10)
        HIGHBD_OBFP_WRAPPER(32, 8, 10)
        HIGHBD_OBFP_WRAPPER(8, 32, 10)
        HIGHBD_OBFP_WRAPPER(16, 4, 10)
        HIGHBD_OBFP_WRAPPER(4, 16, 10)
        break;

      case AOM_BITS_12:
        HIGHBD_BFP_WRAPPER(64, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 64, 12)
        HIGHBD_BFP_WRAPPER(32, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 32, 12)
        HIGHBD_BFP_WRAPPER(16, 4, 12)
        HIGHBD_BFP_WRAPPER(4, 16, 12)
        HIGHBD_BFP_WRAPPER(32, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 32, 12)
        HIGHBD_BFP_WRAPPER(64, 32, 12)
        HIGHBD_BFP_WRAPPER(32, 64, 12)
        HIGHBD_BFP_WRAPPER(32, 32, 12)
        HIGHBD_BFP_WRAPPER(64, 64, 12)
        HIGHBD_BFP_WRAPPER(16, 16, 12)
        HIGHBD_BFP_WRAPPER(16, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 16, 12)
        HIGHBD_BFP_WRAPPER(8, 8, 12)
        HIGHBD_BFP_WRAPPER(8, 4, 12)
        HIGHBD_BFP_WRAPPER(4, 8, 12)
        HIGHBD_BFP_WRAPPER(4, 4, 12)
        HIGHBD_BFP_WRAPPER(128, 128, 12)
        HIGHBD_BFP_WRAPPER(128, 64, 12)
        HIGHBD_BFP_WRAPPER(64, 128, 12)

        HIGHBD_MBFP_WRAPPER(128, 128, 12)
        HIGHBD_MBFP_WRAPPER(128, 64, 12)
        HIGHBD_MBFP_WRAPPER(64, 128, 12)
        HIGHBD_MBFP_WRAPPER(64, 64, 12)
        HIGHBD_MBFP_WRAPPER(64, 32, 12)
        HIGHBD_MBFP_WRAPPER(32, 64, 12)
        HIGHBD_MBFP_WRAPPER(32, 32, 12)
        HIGHBD_MBFP_WRAPPER(32, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 32, 12)
        HIGHBD_MBFP_WRAPPER(16, 16, 12)
        HIGHBD_MBFP_WRAPPER(8, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 8, 12)
        HIGHBD_MBFP_WRAPPER(4, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 4, 12)
        HIGHBD_MBFP_WRAPPER(4, 4, 12)
        HIGHBD_MBFP_WRAPPER(64, 16, 12)
        HIGHBD_MBFP_WRAPPER(16, 64, 12)
        HIGHBD_MBFP_WRAPPER(32, 8, 12)
        HIGHBD_MBFP_WRAPPER(8, 32, 12)
        HIGHBD_MBFP_WRAPPER(16, 4, 12)
        HIGHBD_MBFP_WRAPPER(4, 16, 12)

        HIGHBD_OBFP_WRAPPER(128, 128, 12)
        HIGHBD_OBFP_WRAPPER(128, 64, 12)
        HIGHBD_OBFP_WRAPPER(64, 128, 12)
        HIGHBD_OBFP_WRAPPER(64, 64, 12)
        HIGHBD_OBFP_WRAPPER(64, 32, 12)
        HIGHBD_OBFP_WRAPPER(32, 64, 12)
        HIGHBD_OBFP_WRAPPER(32, 32, 12)
        HIGHBD_OBFP_WRAPPER(32, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 32, 12)
        HIGHBD_OBFP_WRAPPER(16, 16, 12)
        HIGHBD_OBFP_WRAPPER(8, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 8, 12)
        HIGHBD_OBFP_WRAPPER(4, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 4, 12)
        HIGHBD_OBFP_WRAPPER(4, 4, 12)
        HIGHBD_OBFP_WRAPPER(64, 16, 12)
        HIGHBD_OBFP_WRAPPER(16, 64, 12)
        HIGHBD_OBFP_WRAPPER(32, 8, 12)
        HIGHBD_OBFP_WRAPPER(8, 32, 12)
        HIGHBD_OBFP_WRAPPER(16, 4, 12)
        HIGHBD_OBFP_WRAPPER(4, 16, 12)
        break;

      default:
        assert(0 &&
               "cm->seq_params.bit_depth should be AOM_BITS_8, "
               "AOM_BITS_10 or AOM_BITS_12");
    }
  }
}
#endif  // CONFIG_AV1_HIGHBITDEPTH

static AOM_INLINE void variance_partition_alloc(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const int num_64x64_blocks = (cm->seq_params.sb_size == BLOCK_64X64) ? 1 : 4;
  if (cpi->td.vt64x64) {
    if (num_64x64_blocks != cpi->td.num_64x64_blocks) {
      aom_free(cpi->td.vt64x64);
      cpi->td.vt64x64 = NULL;
    }
  }
  if (!cpi->td.vt64x64) {
    CHECK_MEM_ERROR(cm, cpi->td.vt64x64,
                    aom_malloc(sizeof(*cpi->td.vt64x64) * num_64x64_blocks));
    cpi->td.num_64x64_blocks = num_64x64_blocks;
  }
}

static AOM_INLINE void copy_frame_prob_info(AV1_COMP *cpi) {
  FrameProbInfo *const frame_probs = &cpi->frame_probs;
  if (cpi->sf.tx_sf.tx_type_search.prune_tx_type_using_stats) {
    av1_copy(frame_probs->tx_type_probs, default_tx_type_probs);
  }
  if (!cpi->sf.inter_sf.disable_obmc &&
      cpi->sf.inter_sf.prune_obmc_prob_thresh > 0) {
    av1_copy(frame_probs->obmc_probs, default_obmc_probs);
  }
  if (cpi->sf.inter_sf.prune_warped_prob_thresh > 0) {
    av1_copy(frame_probs->warped_probs, default_warped_probs);
  }
  if (cpi->sf.interp_sf.adaptive_interp_filter_search == 2) {
    av1_copy(frame_probs->switchable_interp_probs,
             default_switchable_interp_probs);
  }
}

static AOM_INLINE void restore_cur_buf(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;
  aom_yv12_copy_frame(&cc->copy_buffer, &cm->cur_frame->buf,
                      av1_num_planes(cm));
}

// Coding context that only needs to be restored when recode loop includes
// filtering (deblocking, CDEF, superres post-encode upscale and/or loop
// restoraton).
static AOM_INLINE void restore_extra_coding_context(AV1_COMP *cpi) {
  CODING_CONTEXT *const cc = &cpi->coding_context;
  AV1_COMMON *cm = &cpi->common;
  cm->lf = cc->lf;
  cm->cdef_info = cc->cdef_info;
  cpi->rc = cc->rc;
}

static AOM_INLINE void release_copy_buffer(CODING_CONTEXT *cc) {
  aom_free_frame_buffer(&cc->copy_buffer);
}

static AOM_INLINE int equal_dimensions_and_border(const YV12_BUFFER_CONFIG *a,
                                                  const YV12_BUFFER_CONFIG *b) {
  return a->y_height == b->y_height && a->y_width == b->y_width &&
         a->uv_height == b->uv_height && a->uv_width == b->uv_width &&
         a->y_stride == b->y_stride && a->uv_stride == b->uv_stride &&
         a->border == b->border &&
         (a->flags & YV12_FLAG_HIGHBITDEPTH) ==
             (b->flags & YV12_FLAG_HIGHBITDEPTH);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_ENCODER_UTILS_H_
