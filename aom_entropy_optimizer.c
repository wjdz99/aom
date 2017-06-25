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

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "av1/common/entropymode.h"
#include "av1/common/onyxc_int.h"

#define aom_count_type unsigned int
#define spaces_per_tab 2

FILE *testfile;

static INLINE unsigned int tree_probs_optimizer(unsigned int i,
                                                const aom_tree_index *tree,
                                                const unsigned int *counts,
                                                aom_prob *probs) {
  const int l = tree[i];
  const unsigned int left_count =
      (l <= 0) ? counts[-l] : tree_probs_optimizer(l, tree, counts, probs);
  const int r = tree[i + 1];
  const unsigned int right_count =
      (r <= 0) ? counts[-r] : tree_probs_optimizer(r, tree, counts, probs);
  probs[i >> 1] = get_binary_prob(left_count, right_count);
  return left_count + right_count;
}

static void stats_parser_recursive(aom_count_type **ct_ptr, FILE *probsfile,
                                   int tabs, int dim_of_cts, int *cts_each_dim,
                                   const aom_tree_index *tree,
                                   int flatten_last_dim) {
  if (dim_of_cts < 1) {
    fprintf(stderr, "The dimension of a counts vector should be at least 1!\n");
    exit(EXIT_FAILURE);
  }
  if (dim_of_cts == 1) {
    int k, total_modes = cts_each_dim[0];
    aom_count_type *counts1d = *ct_ptr;
    aom_prob *probs = aom_malloc(sizeof(aom_prob) * (total_modes - 1));

    (*ct_ptr) += total_modes;
    if (tree) {
      tree_probs_optimizer(0, tree, counts1d, probs);
    } else {
      assert(total_modes == 2);
      probs[0] = get_binary_prob(counts1d[0], counts1d[1]);
    }
    if (tabs > 0) fprintf(probsfile, "%*c", tabs * spaces_per_tab, ' ');
    for (k = 0; k < total_modes - 1; ++k) {
      fprintf(probsfile, " %3d,", probs[k]);
      fprintf(testfile, "%d ", counts1d[k]);
    }
    fprintf(testfile, "%d\n", counts1d[k]);
  } else if (dim_of_cts == 2 && flatten_last_dim) {
    int k;

    assert(cts_each_dim[1] == 2);

    for (k = 0; k < cts_each_dim[0]; ++k) {
      fprintf(probsfile, " %3d,", get_binary_prob((*ct_ptr)[0], (*ct_ptr)[1]));
      fprintf(testfile, "%d %d\n", (*ct_ptr)[0], (*ct_ptr)[1]);
      (*ct_ptr) += 2;
    }
  } else {
    int k;

    for (k = 0; k < cts_each_dim[0]; ++k) {
      int tabs_next_level;
      if (dim_of_cts == 2 || (dim_of_cts == 3 && flatten_last_dim)) {
        fprintf(probsfile, "%*c{", tabs * spaces_per_tab, ' ');
        tabs_next_level = 0;
      } else {
        fprintf(probsfile, "%*c{\n", tabs * spaces_per_tab, ' ');
        tabs_next_level = tabs + 1;
      }
      stats_parser_recursive(ct_ptr, probsfile, tabs_next_level, dim_of_cts - 1,
                             cts_each_dim + 1, tree, flatten_last_dim);
      if (dim_of_cts == 2 || (dim_of_cts == 3 && flatten_last_dim))
        fprintf(probsfile, "},\n");
      else
        fprintf(probsfile, "%*c},\n", tabs * spaces_per_tab, ' ');
    }
  }
}

static void stats_parser(aom_count_type *counts, FILE *probsfile,
                         int dim_of_cts, int *cts_each_dim,
                         const aom_tree_index *tree, int flatten_last_dim,
                         char *prefix) {
  aom_count_type *ct_ptr = counts;

  assert(!flatten_last_dim || (cts_each_dim[dim_of_cts - 1] == 2));

  fprintf(probsfile, "%s = {\n", prefix);
  stats_parser_recursive(&ct_ptr, probsfile, 1, dim_of_cts, cts_each_dim, tree,
                         flatten_last_dim);
  fprintf(probsfile, "};\n\n");
  fprintf(testfile, "\n");
}

static void stats_parser_for_ref_frames(FRAME_COUNTS *fc_ptr, FILE *probsfile) {
  int cts_each_dim[10];

  cts_each_dim[0] = COMP_INTER_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc_ptr->comp_inter[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_comp_inter_p"
               "[COMP_INTER_CONTEXTS]");

#if CONFIG_EXT_COMP_REFS
  cts_each_dim[0] = COMP_REF_TYPE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(
      &(fc_ptr->comp_ref_type[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
      "static const aom_prob default_comp_ref_type_p[COMP_REF_TYPE_CONTEXTS]");

  cts_each_dim[0] = UNI_COMP_REF_CONTEXTS;
  cts_each_dim[1] = UNIDIR_COMP_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(
      &(fc_ptr->uni_comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
      "static const aom_prob\n"
      "default_uni_comp_ref_p[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1]");
#endif  // CONFIG_EXT_COMP_REFS

  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = SINGLE_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc_ptr->single_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL,
               1,
               "static const aom_prob default_single_ref_p[REF_CONTEXTS]"
               "[SINGLE_REFS - 1]");

#if CONFIG_EXT_REFS
  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = FWD_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(
      &(fc_ptr->comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
      "static const aom_prob default_comp_ref_p[REF_CONTEXTS][FWD_REFS - 1]");

  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = BWD_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc_ptr->comp_bwdref[0][0][0]), probsfile, 3, cts_each_dim,
               NULL, 1,
               "static const aom_prob "
               "default_comp_bwdref_p[REF_CONTEXTS][BWD_REFS - 1]");
#else
  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = COMP_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc_ptr->comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL,
               1,
               "static const aom_prob default_comp_ref_p[REF_CONTEXTS]"
               "[COMP_REFS - 1]");
#endif  // CONFIG_EXT_REFS
}

int main(int argc, const char **argv_) {
  FILE *statsfile;
  FILE *probsfile = fopen("optimized_probs.c", "w");
  FILE *probsfile_frame_type[FRAME_CONTEXTS];
  int cts_each_dim[10];
  FRAME_COUNTS fc;
  FRAME_COUNTS fc_frame_type[FRAME_CONTEXTS];
  int i;

  testfile = fopen("test.log", "w");
  if (argc < 2) {
    fprintf(stderr, "Please specify the input stats file!\n");
    exit(EXIT_FAILURE);
  }

  statsfile = fopen(argv_[1], "rb");
  if (!statsfile) {
    fprintf(stderr, "Failed to open input file!\n");
    exit(EXIT_FAILURE);
  }

  fread(&fc, sizeof(FRAME_COUNTS), 1, statsfile);
  for (i = 0; i < FRAME_CONTEXTS; ++i) {
    char file_name[256];
    sprintf(file_name, "optimized_probs_type_%d.c", i);
    probsfile_frame_type[i] = fopen(file_name, "w");
    if (!probsfile_frame_type[i]) {
      fprintf(stderr, "Failed to open probs file %s!\n", file_name);
      exit(EXIT_FAILURE);
    }
    fread(&fc_frame_type[i], sizeof(FRAME_COUNTS), 1, statsfile);
  }

  cts_each_dim[0] = INTRA_MODES;
  cts_each_dim[1] = INTRA_MODES;
  cts_each_dim[2] = INTRA_MODES;
  stats_parser(&(fc.kf_y_mode[0][0][0]), probsfile, 3, cts_each_dim,
               av1_intra_mode_tree, 0,
               "const aom_prob av1_kf_y_mode_prob[INTRA_MODES][INTRA_MODES]"
               "[INTRA_MODES - 1]");

  cts_each_dim[0] = BLOCK_SIZE_GROUPS;
  cts_each_dim[1] = INTRA_MODES;
  stats_parser(&(fc.y_mode[0][0]), probsfile, 2, cts_each_dim,
               av1_intra_mode_tree, 0,
               "static const aom_prob default_if_y_probs[BLOCK_SIZE_GROUPS]"
               "[INTRA_MODES - 1]");

  cts_each_dim[0] = INTRA_MODES;
  cts_each_dim[1] = INTRA_MODES;
  stats_parser(&(fc.uv_mode[0][0]), probsfile, 2, cts_each_dim,
               av1_intra_mode_tree, 0,
               "static const aom_prob default_uv_probs[INTRA_MODES]"
               "[INTRA_MODES - 1]");

  cts_each_dim[0] = PARTITION_CONTEXTS;
#if CONFIG_EXT_PARTITION_TYPES
  cts_each_dim[1] = EXT_PARTITION_TYPES;
  // TODO(yuec): Wrong prob for context = 0, because the old tree is used
  stats_parser(&(fc.partition[0][0]), probsfile, 2, cts_each_dim,
               av1_ext_partition_tree, 0,
               "static const aom_prob default_partition_probs"
               "[PARTITION_CONTEXTS][EXT_PARTITION_TYPES - 1]");
#else
  cts_each_dim[1] = PARTITION_TYPES;
  stats_parser(&(fc.partition[0][0]), probsfile, 2, cts_each_dim,
               av1_partition_tree, 0,
               "static const aom_prob default_partition_probs"
               "[PARTITION_CONTEXTS][PARTITION_TYPES - 1]");
#endif

  cts_each_dim[0] = SWITCHABLE_FILTER_CONTEXTS;
  cts_each_dim[1] = SWITCHABLE_FILTERS;
  stats_parser(&(fc.switchable_interp[0][0]), probsfile, 2, cts_each_dim,
               av1_switchable_interp_tree, 0,
               "static const aom_prob \n"
               "default_switchable_interp_prob[SWITCHABLE_FILTER_CONTEXTS]"
               "[SWITCHABLE_FILTERS - 1]");

  cts_each_dim[0] = TX_SIZES;
  cts_each_dim[1] = PLANE_TYPES;
  cts_each_dim[2] = REF_TYPES;
  cts_each_dim[3] = BLOCKZ_CONTEXTS;
  cts_each_dim[4] = 2;
  stats_parser(&(fc.blockz_count[0][0][0][0][0]), probsfile, 5, cts_each_dim,
               NULL, 1,
               "static const aom_prob av1_default_blockzero_probs[TX_SIZES]"
               "[PLANE_TYPES][REF_TYPES][BLOCKZ_CONTEXTS]");

  cts_each_dim[0] = NEWMV_MODE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.newmv_mode[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_newmv_prob[NEWMV_MODE_CONTEXTS]");

  cts_each_dim[0] = ZEROMV_MODE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(
      &(fc.zeromv_mode[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
      "static const aom_prob default_zeromv_prob[ZEROMV_MODE_CONTEXTS]");

  cts_each_dim[0] = REFMV_MODE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.refmv_mode[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_refmv_prob[REFMV_MODE_CONTEXTS]");

  cts_each_dim[0] = DRL_MODE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.drl_mode[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_drl_prob[DRL_MODE_CONTEXTS]");

/*


#if CONFIG_EXT_INTER
  skip_stats(statsfile, INTER_MODE_CONTEXTS * INTER_COMPOUND_MODES);
#if CONFIG_COMPOUND_SINGLEREF
  skip_stats(statsfile, INTER_MODE_CONTEXTS * INTER_SINGLEREF_COMP_MODES);
#endif
#if CONFIG_INTERINTRA
  skip_stats(statsfile, BLOCK_SIZE_GROUPS * 2);
  skip_stats(statsfile, BLOCK_SIZE_GROUPS * INTERINTRA_MODES);
  skip_stats(statsfile, BLOCK_SIZES * 2);
#endif
  skip_stats(statsfile, BLOCK_SIZES * COMPOUND_TYPES);
#endif
*/

#if CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION
  cts_each_dim[0] = BLOCK_SIZES;
  cts_each_dim[1] = MOTION_MODES;
  stats_parser(&(fc.motion_mode[0][0]), probsfile, 2, cts_each_dim,
               av1_motion_mode_tree, 0,
               "static const aom_prob default_motion_mode_prob[BLOCK_SIZES]"
               "[MOTION_MODES - 1]");
#if CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
  cts_each_dim[0] = BLOCK_SIZES;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.obmc[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_obmc_prob[BLOCK_SIZES]");
#endif  // CONFIG_MOTION_VAR && CONFIG_WARPED_MOTION
#endif  // CONFIG_MOTION_VAR || CONFIG_WARPED_MOTION

  cts_each_dim[0] = INTRA_INTER_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.intra_inter[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_intra_inter_p"
               "[INTRA_INTER_CONTEXTS]");

  stats_parser_for_ref_frames(&fc, probsfile);
  for (i = 0; i < FRAME_CONTEXTS; ++i)
    stats_parser_for_ref_frames(&fc_frame_type[i], probsfile_frame_type[i]);

  /*
  #if CONFIG_EXT_INTER && CONFIG_COMPOUND_SINGLEREF
    skip_stats(statsfile, COMP_INTER_MODE_CONTEXTS * 2);
  #endif

    // TODO(yuec): move tx_size_totals to where only encoder will use
    skip_stats(statsfile, TX_SIZES);
    // TODO(yuec): av1_tx_size_tree has variable size
    skip_stats(statsfile, MAX_TX_DEPTH * TX_SIZE_CONTEXTS * (MAX_TX_DEPTH + 1));
  #if CONFIG_EXT_TX && CONFIG_RECT_TX && CONFIG_RECT_TX_EXT
    skip_stats(statsfile, 2);
  #endif
  #if CONFIG_VAR_TX
    skip_stats(statsfile, TXFM_PARTITION_CONTEXTS * 2);
  #endif
  */

  cts_each_dim[0] = SKIP_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.skip[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_skip_probs[SKIP_CONTEXTS]");
/*

  skip_stats(statsfile, (sizeof(nmv_context_counts) / sizeof(aom_count_type)) *
                            NMV_CONTEXTS);

#if CONFIG_INTRABC
  skip_stats(statsfile, 2);
  skip_stats(statsfile, sizeof(nmv_context_counts) / sizeof(aom_count_type));
#endif

#if CONFIG_DELTA_Q
  skip_stats(statsfile, DELTA_Q_PROBS * 2);
#if CONFIG_EXT_DELTA_Q
  skip_stats(statsfile, DELTA_LF_PROBS * 2);
#endif
#endif
*/

#if CONFIG_EXT_TX
/*
#if CONFIG_RECT_TX
skip_stats(statsfile, TX_SIZES * TX_SIZES);
#endif
skip_stats(statsfile, EXT_TX_SETS_INTER * EXT_TX_SIZES * TX_TYPES);
skip_stats(statsfile,
           EXT_TX_SETS_INTRA * EXT_TX_SIZES * INTRA_MODES * TX_TYPES);
           */
#else
  /*  // TODO(yuec): intra_ext_tx use different trees depending on the context
    skip_stats(statsfile, EXT_TX_SIZES * TX_TYPES * TX_TYPES);*/

  cts_each_dim[0] = EXT_TX_SIZES;
  cts_each_dim[1] = TX_TYPES;
  stats_parser(&(fc.inter_ext_tx[0][0]), probsfile, 2, cts_each_dim,
               av1_ext_tx_tree, 0,
               "static const aom_prob default_inter_ext_tx_prob"
               "[EXT_TX_SIZES][TX_TYPES - 1]");
#endif
  /*

  #if CONFIG_SUPERTX
    skip_stats(statsfile, PARTITION_SUPERTX_CONTEXTS * TX_SIZES * 2);
    skip_stats(statsfile, TX_SIZES);
  #endif

    skip_stats(statsfile, sizeof(struct seg_counts) / sizeof(aom_count_type));

  #if CONFIG_EXT_INTRA
  #if CONFIG_INTRA_INTERP
    skip_stats(statsfile, (INTRA_FILTERS + 1) * INTRA_FILTERS);
  #endif
  #endif

  #if CONFIG_FILTER_INTRA
    skip_stats(statsfile, PLANE_TYPES * 2);
  #endif
  */

  fclose(statsfile);
  fclose(testfile);
  fclose(probsfile);
  for (i = 0; i < FRAME_CONTEXTS; ++i) fclose(probsfile_frame_type[i]);

  return 1;
}
