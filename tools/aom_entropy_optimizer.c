/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <assert.h>
#include <stdio.h>
#include "./aom_config.h"
#include "av1/common/entropymode.h"

#define spaces_per_tab 2

typedef unsigned int aom_count_type;
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

    if (probs == NULL) {
      fprintf(stderr, "Allocating prob array fails!\n");
      exit(EXIT_FAILURE);
    }

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

int main(int argc, const char **argv_) {
  FILE *statsfile;
  FILE *probsfile = fopen("optimized_probs.c", "w");
  int cts_each_dim[10];
  FRAME_COUNTS fc;

  testfile = fopen("aom_entropy_optimizer_parsed_counts.log", "w");
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

#if CONFIG_EXT_INTER
  cts_each_dim[0] = INTER_MODE_CONTEXTS;
  cts_each_dim[1] = INTER_COMPOUND_MODES;
  stats_parser(&(fc.inter_compound_mode[0][0]), probsfile, 2, cts_each_dim,
               av1_inter_compound_mode_tree, 0,
               "static const aom_prob default_inter_compound_mode_probs\n"
               "[INTER_MODE_CONTEXTS][INTER_COMPOUND_MODES - 1]");
#if CONFIG_COMPOUND_SINGLEREF
  cts_each_dim[0] = INTER_MODE_CONTEXTS;
  cts_each_dim[1] = INTER_SINGLEREF_COMP_MODES;
  stats_parser(&(fc.inter_singleref_comp_mode[0][0]), probsfile, 2,
               cts_each_dim, av1_inter_singleref_comp_mode_tree, 0,
               "static const aom_prob default_inter_singleref_comp_mode_probs\n"
               "[INTER_MODE_CONTEXTS][INTER_SINGLEREF_COMP_MODES - 1]");
#endif
#if CONFIG_INTERINTRA
  cts_each_dim[0] = BLOCK_SIZE_GROUPS;
  cts_each_dim[1] = 2;
  stats_parser(
      &(fc.interintra[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
      "static const aom_prob default_interintra_prob[BLOCK_SIZE_GROUPS]");

  cts_each_dim[0] = BLOCK_SIZE_GROUPS;
  cts_each_dim[1] = INTERINTRA_MODES;
  stats_parser(
      &(fc.interintra_mode[0][0]), probsfile, 2, cts_each_dim,
      av1_interintra_mode_tree, 0,
      "static const aom_prob"
      "default_interintra_mode_prob[BLOCK_SIZE_GROUPS][INTERINTRA_MODES - 1]");

  cts_each_dim[0] = BLOCK_SIZES;
  cts_each_dim[1] = 2;
  stats_parser(
      &(fc.wedge_interintra[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
      "static const aom_prob default_wedge_interintra_prob[BLOCK_SIZES]");
#endif
  cts_each_dim[0] = BLOCK_SIZES;
  cts_each_dim[1] = COMPOUND_TYPES;
  stats_parser(&(fc.compound_interinter[0][0]), probsfile, 2, cts_each_dim,
               av1_compound_type_tree, 0,
               "static const aom_prob default_compound_type_probs"
               "[BLOCK_SIZES]TYPES - 1]");
#endif

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

  cts_each_dim[0] = COMP_INTER_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.comp_inter[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_comp_inter_p"
               "[COMP_INTER_CONTEXTS]");

#if CONFIG_EXT_COMP_REFS
  cts_each_dim[0] = COMP_REF_TYPE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(
      &(fc.comp_ref_type[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
      "static const aom_prob default_comp_ref_type_p[COMP_REF_TYPE_CONTEXTS]");

  cts_each_dim[0] = UNI_COMP_REF_CONTEXTS;
  cts_each_dim[1] = UNIDIR_COMP_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(
      &(fc.uni_comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
      "static const aom_prob\n"
      "default_uni_comp_ref_p[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS - 1]");
#endif

  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = SINGLE_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc.single_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
               "static const aom_prob default_single_ref_p[REF_CONTEXTS]"
               "[SINGLE_REFS - 1]");

#if CONFIG_EXT_REFS
  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = FWD_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(
      &(fc.comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
      "static const aom_prob default_comp_ref_p[REF_CONTEXTS][FWD_REFS - 1]");

  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = BWD_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc.comp_bwdref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
               "static const aom_prob "
               "default_comp_bwdref_p[REF_CONTEXTS][BWD_REFS - 1]");
#else
  cts_each_dim[0] = REF_CONTEXTS;
  cts_each_dim[1] = COMP_REFS - 1;
  cts_each_dim[2] = 2;
  stats_parser(&(fc.comp_ref[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
               "static const aom_prob default_comp_ref_p[REF_CONTEXTS]"
               "[COMP_REFS - 1]");
#endif  // CONFIG_EXT_REFS

#if CONFIG_EXT_INTER && CONFIG_COMPOUND_SINGLEREF
  cts_each_dim[0] = COMP_INTER_MODE_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.comp_inter_mode[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob "
               "default_comp_inter_mode_p[COMP_INTER_MODE_CONTEXTS]");
#endif

// TODO(yuec): av1_tx_size_tree has variable sizes, so needs special handling
#if CONFIG_EXT_TX && CONFIG_RECT_TX && CONFIG_RECT_TX_EXT
  cts_each_dim[0] = 2;
  stats_parser(&(fc.quarter_tx_size[0]), probsfile, 1, cts_each_dim, NULL, 1,
               "static const aom_prob default_quarter_tx_size_prob");
#endif
#if CONFIG_VAR_TX
  cts_each_dim[0] = TXFM_PARTITION_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.txfm_partition[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob "
               "default_txfm_partition_probs[TXFM_PARTITION_CONTEXTS]");
#endif

  cts_each_dim[0] = SKIP_CONTEXTS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.skip[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_skip_probs[SKIP_CONTEXTS]");

#if CONFIG_INTRABC
  cts_each_dim[0] = 2;
  stats_parser(&(fc.intrabc[0]), probsfile, 1, cts_each_dim, NULL, 1,
               "INTRABC_PROB_DEFAULT");
#endif

#if CONFIG_DELTA_Q
  cts_each_dim[0] = DELTA_Q_PROBS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.delta_q[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_delta_q_probs[DELTA_Q_PROBS]");
#if CONFIG_EXT_DELTA_Q
  cts_each_dim[0] = DELTA_LF_PROBS;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.delta_lf[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_delta_lf_probs[DELTA_LF_PROBS]");
#endif
#endif

#if CONFIG_EXT_TX
// TODO(yuec): different trees are used depending on selected ext tx set
#else
  // TODO(yuec): intra_ext_tx use different trees depending on the context
  cts_each_dim[0] = EXT_TX_SIZES;
  cts_each_dim[1] = TX_TYPES;
  stats_parser(&(fc.inter_ext_tx[0][0]), probsfile, 2, cts_each_dim,
               av1_ext_tx_tree, 0,
               "static const aom_prob default_inter_ext_tx_prob"
               "[EXT_TX_SIZES][TX_TYPES - 1]");
#endif

#if CONFIG_SUPERTX
  cts_each_dim[0] = PARTITION_SUPERTX_CONTEXTS;
  cts_each_dim[1] = TX_SIZES;
  cts_each_dim[2] = 2;
  stats_parser(&(fc.supertx[0][0][0]), probsfile, 3, cts_each_dim, NULL, 1,
               "static const aom_prob\n"
               "default_supertx_prob[PARTITION_SUPERTX_CONTEXTS][TX_SIZES]");
#endif

#if CONFIG_EXT_INTRA
#if CONFIG_INTRA_INTERP
  cts_each_dim[0] = INTRA_FILTERS + 1;
  cts_each_dim[1] = INTRA_FILTERS;
  stats_parser(
      &(fc.intra_filter[0][0]), probsfile, 2, cts_each_dim,
      av1_intra_filter_tree, 0,
      "static const aom_prob\n"
      "default_intra_filter_probs[INTRA_FILTERS + 1][INTRA_FILTERS - 1]");
#endif
#endif

#if CONFIG_FILTER_INTRA
  cts_each_dim[0] = PLANE_TYPES;
  cts_each_dim[1] = 2;
  stats_parser(&(fc.filter_intra[0][0]), probsfile, 2, cts_each_dim, NULL, 1,
               "static const aom_prob default_filter_intra_probs[2]");
#endif

  fclose(statsfile);
  fclose(testfile);
  fclose(probsfile);

  return 1;
}
