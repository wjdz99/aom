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
#ifndef AOM_AV1_ENCODER_PICKRST_H_
#define AOM_AV1_ENCODER_PICKRST_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/encoder.h"
#include "aom_ports/system_state.h"

struct yv12_buffer_config;
struct AV1_COMP;

static const uint8_t g_shuffle_stats_data[16] = {
  0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8,
};

static const uint8_t g_shuffle_stats_highbd_data[32] = {
  0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7, 8, 9,
  0, 1, 2, 3, 2, 3, 4, 5, 4, 5, 6, 7, 6, 7, 8, 9,
};

static INLINE uint8_t find_average(const uint8_t *src, int h_start, int h_end,
                                   int v_start, int v_end, int stride) {
  uint64_t sum = 0;
  for (int i = v_start; i < v_end; i++) {
    for (int j = h_start; j < h_end; j++) {
      sum += src[i * stride + j];
    }
  }
  uint64_t avg = sum / ((v_end - v_start) * (h_end - h_start));
  return (uint8_t)avg;
}

static INLINE uint16_t find_average_highbd(const uint16_t *src, int h_start,
                                           int h_end, int v_start, int v_end,
                                           int stride) {
  uint64_t sum = 0;
  for (int i = v_start; i < v_end; i++) {
    for (int j = h_start; j < h_end; j++) {
      sum += src[i * stride + j];
    }
  }
  uint64_t avg = sum / ((v_end - v_start) * (h_end - h_start));
  return (uint16_t)avg;
}

void av1_pick_filter_restoration(const YV12_BUFFER_CONFIG *sd,
#if CONFIG_LOOP_RESTORE_CNN
                                 bool allow_restore_cnn_y,
#endif  // CONFIG_LOOP_RESTORE_CNN
                                 AV1_COMP *cpi);

// Function type to determine edge cost - trivial if cost is constant, but
// necessary when edge cost is dependent on path from src
// info : pointer to unspecified structure, holds any information needed
//  to calculate edge cost
// path :  pointer to Vector holding current path to edge
// node_idx : start of edge
// max_out_nodes: max outgoing edges. Without subsets, this is the number of
//  nodes in the graph. Otherwise, this is the number of nodes in each subset.
// int out_edge: outgoing edge we are calculating cost for, numbered in
//  relation to node_idx
typedef double (*graph_edge_cost_t)(void *info, Vector *path, int node_idx,
                                    int max_out_nodes, int out_edge);

// Searching a directed graph.
// If subsets == false, denotes a graph where nodes can have outgoing edge to
// any other node.
// Otherwise, nodes are organized into equally sized subsets. Any outgoing
// edges from nodes in a subset only go to nodes in one other subset, and
// there are no cycles between subsets. We are finding a set of nodes, one per
// subset, that form a min cost path from src to dest.
// (Use case is finding restoration types for each unit in RESTORE_SWITCHABLE.)
// node_idx : start of currently explored path
// dest_idx : destination of currently explored path
// max_out_nodes: max outgoing edges. Without subsets, this is the number of
//  nodes in the graph. Otherwise, this is the number of nodes in each subset.
// graph: pointer to adjacency matrix to indicate edges between nodes. If no
//  edge is present between nodes, element is set to INFINITY. Without subsets,
//  this is an n*n matrix where n is the number of nodes and element
//  graph[n1][n2] is an edge from n1 to n2. With subsets, this is a
//  (# of subsets * max_out_nodes + 2 nodes for src and dest) * max_out_nodes
//  matrix and element graph [n][e] is an edge from n to the eth node in the
//  next subset (see rst_mergecoeffs_tst.c for example)
// best_path : pointer to Vector storing best path from start to destination
// subsets : indicates whether graph needs to be organized into subsets
// cost_fn : function to dynamically determine edge cost
// info : pointer to unspecified structure, holds any information needed
//  to calculate edge cost
double min_cost_type_path(int node_idx, int dest_idx, int max_out_nodes,
                          const double *graph, Vector *best_path, bool subsets,
                          graph_edge_cost_t cost_fn, void *info);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_PICKRST_H_
