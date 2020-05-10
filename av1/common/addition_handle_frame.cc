/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/common/addition_handle_frame.h"

uint8_t *getYbuf(uint8_t *yPxl, int height, int width, int stride) {
  if (!yPxl) return nullptr;
  int buflen = height * width;
  unsigned char *buf = new unsigned char[buflen];  //定义一块一帧图像大小的buf
  unsigned char *p = buf;
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      unsigned char uctemp = (unsigned char)(*(yPxl + x));
      *p = uctemp;
      p++;
    }
    yPxl += stride;  //一行的一个Y值对应下一行同样位置的增量
  }
  return buf;
}

/*Feed full frame image into the network*/
void addition_handle_frame(AV1_COMMON *cm) {
  YV12_BUFFER_CONFIG *pcPicYuvRec = &cm->cur_frame->buf;

  if (!cm->seq_params.use_highbitdepth) {
    uint8_t *py = pcPicYuvRec->y_buffer;
    uint8_t *bkuPy = py;

    int height = pcPicYuvRec->y_height;
    int width = pcPicYuvRec->y_width;
    int stride = pcPicYuvRec->y_stride;
    uint8_t **buf = new uint8_t *[height];
    for (int i = 0; i < height; i++) {
      buf[i] = new uint8_t[width];
    }
    // if (cm->current_frame.order_hint % 16 == 0 &&
    //    cm->current_frame.order_hint != 0)
    //  buf = TF_Predict(py, height, width, stride, cm->base_qindex - 80,
    //                   cm->current_frame.frame_type);
    // else
    buf = TF_Predict(py, height, width, stride, cm->base_qindex,
                     cm->current_frame.frame_type);
    /* buf = TF_Predict(py, height, width, stride,
     * cm->current_frame.frame_type==0?170:212);*/
    // buf = TF_Predict(py, height, width, stride, 2);

    // FILE *ff = fopen("D:/AOMedia/aom_org/test_result/cnn1_cnn2.yuv", "wb");
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        *(bkuPy + j) = buf[i][j];  // Fill in the luma buffer again
                                   // fwrite(bkuPy + j, 1, 1, ff);
      }
      bkuPy += stride;
    }
    // fclose(ff);
  } else {
    uint16_t *py = CONVERT_TO_SHORTPTR(pcPicYuvRec->y_buffer);
    uint16_t *bkuPy = py;

    int height = pcPicYuvRec->y_height;
    int width = pcPicYuvRec->y_width;
    int stride = pcPicYuvRec->y_stride;

    uint16_t **buf = new uint16_t *[height];
    for (int i = 0; i < height; i++) {
      buf[i] = new uint16_t[width];
    }

    buf = TF_Predict_hbd(py, height, width, stride);

    // FILE *ff = fopen("end.yuv", "wb");
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        *(bkuPy + j) = buf[i][j];  // Fill in the luma buffer again
        // fwrite(bkuPy + j, 1, 1, ff);
      }
      bkuPy += stride;
    }
  }
  finish_python();
}


/*Split into 1000x1000 blocks into the network*/
void addition_handle_blocks(AV1_COMMON *cm) {
  YV12_BUFFER_CONFIG *pcPicYuvRec = &cm->cur_frame->buf;

  int height = pcPicYuvRec->y_height;
  int width = pcPicYuvRec->y_width;
  int buf_width = 1000;
  int buf_height = 1000;
  int bufx_count = width / buf_width;
  int bufy_count = height / buf_height;
  int stride = pcPicYuvRec->y_stride;
  int cur_buf_width;
  int cur_buf_height;

  init_python();

  if (!cm->seq_params.use_highbitdepth) {
    uint8_t *py = pcPicYuvRec->y_buffer;
    uint8_t *bkuPy = py;
    uint8_t *bkuPyTemp = bkuPy;

    for (int y = 0; y < bufy_count + 1; y++) {
      for (int x = 0; x < bufx_count + 1; x++) {
        if (x == bufx_count)
          cur_buf_width = width - buf_width * bufx_count;
        else
          cur_buf_width = buf_width;
        if (y == bufy_count)
          cur_buf_height = height - buf_height * bufy_count;
        else
          cur_buf_height = buf_height;

        if (cur_buf_width != 0 && cur_buf_height != 0) {
          uint8_t **buf = new uint8_t *[cur_buf_height];
          for (int i = 0; i < cur_buf_height; i++) {
            buf[i] = new uint8_t[cur_buf_width];
          }

          TF_Predict_block(buf, py + x * buf_width + stride * buf_height * y,
                           cur_buf_height, cur_buf_width, stride,
                           cm->base_qindex);

          bkuPy = bkuPyTemp;
          // FILE *ff = fopen("end.yuv", "wb");
          for (int i = 0; i < cur_buf_height; i++) {
            for (int j = 0; j < cur_buf_width; j++) {
              *(bkuPy + j + x * buf_width + y * buf_height * stride) =
                  buf[i][j];  // Fill in the luma buffer again
              // fwrite(bkuPy + j, 1, 1, ff);
            }
            bkuPy += stride;
          }
          // fclose(ff);
        }
      }
    }
  } else {
    uint16_t *py = CONVERT_TO_SHORTPTR(pcPicYuvRec->y_buffer);
    uint16_t *bkuPy = py;
    uint16_t *bkuPyTemp = bkuPy;

    for (int y = 0; y < bufy_count + 1; y++) {
      for (int x = 0; x < bufx_count + 1; x++) {
        if (x == bufx_count)
          cur_buf_width = width - buf_width * bufx_count;
        else
          cur_buf_width = buf_width;
        if (y == bufy_count)
          cur_buf_height = height - buf_height * bufy_count;
        else
          cur_buf_height = buf_height;

        if (cur_buf_width != 0 && cur_buf_height != 0) {
          uint16_t **buf = new uint16_t *[cur_buf_height];
          for (int i = 0; i < cur_buf_height; i++) {
            buf[i] = new uint16_t[cur_buf_width];
          }

          buf =
              TF_Predict_block_hbd(py + x * buf_width + stride * buf_height * y,
                                   cur_buf_height, cur_buf_width, stride);

          bkuPy = bkuPyTemp;
          // FILE *ff = fopen("end.yuv", "wb");
          for (int i = 0; i < cur_buf_height; i++) {
            for (int j = 0; j < cur_buf_width; j++) {
              *(bkuPy + j + x * buf_width + y * buf_height * stride) =
                  buf[i][j];
              // fwrite(bkuPy + j, 1, 1, ff);
            }
            bkuPy += stride;
          }
          // fclose(ff);
        }
      }
    }
  }
  finish_python();
}

#if SAVE_PRED
void addition_handle_blocks_adp(AV1_COMMON *cm, double threshold_MAD,
                                int block_size) {
#define out_yuv 1
  YV12_BUFFER_CONFIG *recon_buf = &cm->cur_frame->buf;
  YV12_BUFFER_CONFIG *pred_buf = &cm->cur_frame->pred_buf;
  int height = recon_buf->y_height;
  int width = recon_buf->y_width;
  int stride = recon_buf->y_stride;

  /*char unfilter_file[256] =
      "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
      "result\\E19080601/unfilter.yuv";
  char filter_file[256] =
      "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
      "result\\E19080601/filter.yuv";
  save_frame_yuv(unfilter_file, recon_buf, cm->current_frame.order_hint);*/
  if (height != 0 && width != 0) {
    int16_t **y_resi_buf = new int16_t *[height];
#if out_yuv
    uint8_t **y_recon_buf = new uint8_t *[height];
#endif
    for (int i = 0; i < height; i++) {
      y_resi_buf[i] = new int16_t[width];
#if out_yuv
      y_recon_buf[i] = new uint8_t[width];
#endif
      for (int j = 0; j < width; j++) {
        y_resi_buf[i][j] = *(recon_buf->y_buffer + j + i * stride) -
                           *(pred_buf->y_buffer + j + i * stride);
        // printf("\ nresi[%d][%d]:%d",i,j, y_resi_buf[i][j]);
      }
    }
#if out_yuv
    uint8_t **u_recon_buf = new uint8_t *[height / 2];
    uint8_t **v_recon_buf = new uint8_t *[height / 2];
    for (int i = 0; i < height / 2; i++) {
      u_recon_buf[i] = new uint8_t[width / 2];
      v_recon_buf[i] = new uint8_t[width / 2];
    }
#endif
    int buf_width = block_size;
    int buf_height = block_size;
    int bufx_count = width / buf_width;
    int bufy_count = height / buf_height;

    int padding_size = 8;

    int cur_buf_width;
    int cur_buf_height;

    int cur_block_width;
    int cur_block_height;

    int start_read_x = 0, start_read_y = 0;
    int start_write_x = 0, start_write_y = 0;

    // if (!cm->seq_params.use_highbitdepth) {
    uint8_t *py = recon_buf->y_buffer;
    uint8_t *bkuPy = py;
    uint8_t *bkuPyTemp = bkuPy;
    for (int y = 0; y < bufy_count + 1; y++) {
      for (int x = 0; x < bufx_count + 1; x++) {
        if (x == bufx_count) {
          cur_block_width = width - buf_width * bufx_count;
          cur_buf_width = width - buf_width * bufx_count + padding_size;
          start_read_x = buf_width * bufx_count - padding_size;
          start_write_x = padding_size;
        } else if (x * buf_width < padding_size) {
          cur_block_width = buf_width;
          cur_buf_width = buf_width + padding_size;
          start_read_x = x * buf_width;
          start_write_x = x * buf_width;
        } else {
          cur_block_width = buf_width;
          cur_buf_width = buf_width + 2 * padding_size;
          start_read_x = x * buf_width - padding_size;
          start_write_x = padding_size;
        }

        if (y == bufy_count) {
          cur_block_height = height - buf_height * bufy_count;
          cur_buf_height = height - buf_height * bufy_count + padding_size;
          start_read_y = buf_height * bufy_count - padding_size;
          start_write_y = padding_size;
        } else if (y * buf_height < padding_size) {
          cur_block_height = buf_height;
          cur_buf_height = buf_height + padding_size;
          start_read_y = y * buf_height;
          start_write_y = y * buf_height;
        } else {
          cur_block_height = buf_height;
          cur_buf_height = buf_height + 2 * padding_size;
          start_read_y = y * buf_height - padding_size;
          start_write_y = padding_size;
        }

        if (cur_buf_width != 0 && cur_buf_height != 0) {
          uint8_t **buf = new uint8_t *[cur_buf_height];
          for (int i = 0; i < cur_buf_height; i++) {
            buf[i] = new uint8_t[cur_buf_width];
          }
          double cur_block_SAD = 0, cur_block_MAD = 0;
          for (int i = 0; i < cur_block_height; i++) {
            for (int j = 0; j < cur_block_width; j++) {
              /* printf("\ nresi[%d][%d]:%d", i + y * buf_height,
                      j + x * buf_width,
                      y_resi_buf[i + y * buf_height][j + x * buf_width]);*/
              cur_block_SAD +=
                  abs(y_resi_buf[i + y * buf_height][j + x * buf_width]);
            }
          }
          cur_block_MAD = cur_block_SAD / cur_block_height / cur_block_width;
          /*  printf("\n cur_block_SAD = %lf, cur_block_MAD = %lf",
             cur_block_SAD, cur_block_MAD);*/

          uint8_t **filter_block_buf = new uint8_t *[cur_block_height];
          for (int i = 0; i < cur_block_height; i++) {
            filter_block_buf[i] = new uint8_t[cur_block_width];
          }
          if (cur_block_MAD > threshold_MAD) {
            // printf("\n x:%d,y:%d", x, y);
            if (cm->current_frame.order_hint % 16 == 0 &&
                cm->current_frame.order_hint != 0 && cm->base_qindex == 212)
              buf = TF_Predict(py + start_read_x + stride * start_read_y,
                               cur_buf_height, cur_buf_width, stride,
                               cm->base_qindex - 80,
                               cm->current_frame.frame_type);
            else {
              buf = TF_Predict(py + start_read_x + stride * start_read_y,
                               cur_buf_height, cur_buf_width, stride,
                               cm->base_qindex, cm->current_frame.frame_type);
            }
            for (int i = 0; i < cur_block_height; i++) {
              for (int j = 0; j < cur_block_width; j++) {
                *(recon_buf->y_buffer + (i + y * buf_height) * stride + j +
                  x * buf_width) = buf[i + start_write_y][j + start_write_x];
#if out_yuv
                y_recon_buf[y * buf_height + i][x * buf_width + j] =
                    buf[i + start_write_y][j + start_write_x];
                u_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                    13;
                v_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                    165;
#endif
              }
            }
          }
#if out_yuv
          else {
            for (int i = 0; i < cur_block_height; i++) {
              for (int j = 0; j < cur_block_width; j++) {
                y_recon_buf[y * buf_height + i][x * buf_width + j] =
                    *(recon_buf->y_buffer + (i + y * buf_height) * stride + j +
                      x * buf_width);
                u_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                    *(recon_buf->u_buffer +
                      (i + y * buf_height) / 2 * recon_buf->uv_stride +
                      (j + x * buf_width) / 2);
                v_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                    *(recon_buf->v_buffer +
                      (i + y * buf_height) / 2 * recon_buf->uv_stride +
                      (j + x * buf_width) / 2);
              }
            }
          }
#endif
        }
      }
    }
#if out_yuv
    char recon_file[256];
    sprintf(recon_file,
            "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
            "result\\E19081201"
            "\\BlowingBubbles_416x240_32x32_adp.yuv");

    save_buf_to_yuv(y_recon_buf, height, width, recon_file,
                    cm->cur_frame->order_hint);
    save_buf_to_yuv(u_recon_buf, height / 2, width / 2, recon_file, 1);
    save_buf_to_yuv(v_recon_buf, height / 2, width / 2, recon_file, 1);
#endif
    // for (int i = 0; i < height; i++) {
    //  for (int j = 0; j < width; j++) {
    //    *(bkuPy + j) = y_filter_buf[i][j];  // fill in the luma buffer again
    //  }

    //  bkuPy += stride;
    //}
  }
  // save_frame_yuv(filter_file, recon_buf, cm->current_frame.order_hint);
}
#endif

// 20191231 CRLC block
#if CRLC_LF
void addition_handle_blocks_CRLC(YV12_BUFFER_CONFIG *source_frame,
                                 AV1_COMMON *cm) {
#define out_yuv 0
  YV12_BUFFER_CONFIG *recon_buf = &cm->cur_frame->buf;
  YV12_BUFFER_CONFIG *source_buf = source_frame;
  int height = recon_buf->y_height;
  int width = recon_buf->y_width;
  int dgr_stride = recon_buf->y_stride;
  uint8_t *y_source_buf = getYbuf(source_buf->y_buffer, source_buf->y_height,
                                  source_buf->y_width, source_buf->y_stride);
  // printf("%d,", y_source_buf[0]);
  uint8_t *y_unfilter_buf =
      getYbuf(recon_buf->y_buffer, height, width, dgr_stride);
  /* char unfilter_file[256] =
       "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
       "result\\E19080601/unfilter.yuv";
   char filter_file[256] =
       "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
       "result\\E19080601/filter.yuv";
   save_frame_yuv(unfilter_file, recon_buf, cm->current_frame.order_hint);*/
  if (height != 0 && width != 0) {
    uint8_t **y_filter_buf = new uint8_t *[height];
#if out_yuv
    uint8_t **y_recon_buf = new uint8_t *[height];
#endif
    for (int i = 0; i < height; i++) {
      y_filter_buf[i] = new uint8_t[width];
#if out_yuv
      y_recon_buf[i] = new uint8_t[width];
#endif
    }
#if out_yuv
    uint8_t **u_recon_buf = new uint8_t *[height / 2];
    uint8_t **v_recon_buf = new uint8_t *[height / 2];
    for (int i = 0; i < height / 2; i++) {
      u_recon_buf[i] = new uint8_t[width / 2];
      v_recon_buf[i] = new uint8_t[width / 2];
    }
#endif
    int buf_width = cm->rst_info[0].restoration_unit_size;
    int buf_height = cm->rst_info[0].restoration_unit_size;
    int bufx_count = AOMMAX((width - 1) / buf_width + 1, 1);

    int bufy_count = AOMMAX((height - 1) / buf_height + 1, 1);

    int padding_size = 8;

    int cur_buf_width = 0;
    int cur_buf_height = 0;

    int cur_block_width = 0;
    int cur_block_height = 0;

    int start_read_x = 0, start_read_y = 0;
    int start_write_x = 0, start_write_y = 0;

    uint8_t *dgrPtr = recon_buf->y_buffer;
    uint8_t *bkuPy = dgrPtr;
    uint8_t *srcPtr = source_buf->y_buffer;
    for (int y = 0; y < bufy_count; y++) {
      for (int x = 0; x < bufx_count; x++) {
        if (x == bufx_count - 1) {
          cur_block_width = width - buf_width * (bufx_count - 1);
          cur_buf_width = width - buf_width * (bufx_count - 1) + padding_size;
          start_read_x = buf_width * (bufx_count - 1) - padding_size;
          start_write_x = padding_size;
        } else if (x * buf_width < padding_size) {
          cur_block_width = buf_width;
          cur_buf_width = buf_width + padding_size;
          start_read_x = x * buf_width;
          start_write_x = x * buf_width;
        } else {
          cur_block_width = buf_width;
          cur_buf_width = buf_width + 2 * padding_size;
          start_read_x = x * buf_width - padding_size;
          start_write_x = padding_size;
        }

        if (y == bufy_count - 1) {
          cur_block_height = height - buf_height * (bufy_count - 1);
          cur_buf_height =
              height - buf_height * (bufy_count - 1) + padding_size;
          start_read_y = buf_height * (bufy_count - 1) - padding_size;
          start_write_y = padding_size;
        } else if (y * buf_height < padding_size) {
          cur_block_height = buf_height;
          cur_buf_height = buf_height + padding_size;
          start_read_y = y * buf_height;
          start_write_y = y * buf_height;
        } else {
          cur_block_height = buf_height;
          cur_buf_height = buf_height + 2 * padding_size;
          start_read_y = y * buf_height - padding_size;
          start_write_y = padding_size;
        }

        if (cur_buf_width != 0 && cur_buf_height != 0) {
          uint8_t **buf = new uint8_t *[cur_buf_height];
          for (int i = 0; i < cur_buf_height; i++) {
            buf[i] = new uint8_t[cur_buf_width];
          }

          uint8_t **filter_block_buf = new uint8_t *[cur_block_height];
          for (int i = 0; i < cur_block_height; i++) {
            filter_block_buf[i] = new uint8_t[cur_block_width];
          }
          /*uint8_t *A[2];*/

          buf =
              TF_Predict_CRLC(dgrPtr + start_read_x + dgr_stride * start_read_y,
                              srcPtr + start_read_x + dgr_stride * start_read_y,
                              &cm->crlc_info[0], cur_buf_height, cur_buf_width,
                              dgr_stride, source_frame->y_stride,
                              cm->base_qindex, cm->current_frame.frame_type);
          for (int i = 0; i < cur_block_height; i++) {
            for (int j = 0; j < cur_block_width; j++) {
              filter_block_buf[i][j] =
                  buf[i + start_write_y][j + start_write_x];
            }
          }
          //printf(
          //    "idx:%d , A0:%u\n", y * bufx_count + x,
          //    cm->rst_info[0].unit_info[y * bufx_count + x].crlc_info.xqd[0]);
          // printf(
          //    "CRLC1:%u\n",
          //    cm->rst_info[0].unit_info[x * bufy_count + x].crlc_info.xqd[1]);

          /*  char source_block_file[128];
          sprintf(source_block_file,
                  "D:\\AOMedia\\aom_org\\test_result\\E19072501\\block_"
                  "yuv\\source_block_%d_%d_%dx%d.yuv",
                  x, y, cur_block_width, cur_block_height);
          save_buf_to_yuv(source_block_buf, cur_block_height, cur_block_width,
                          source_block_file);
          char unfilter_block_file[128];
          sprintf(unfilter_block_file,
                  "D:\\AOMedia\\aom_org\\test_result\\E19072501\\block_"
                  "yuv\\unfilter_block_%"
                  "d_%d_%dx%d.yuv",
                  x, y, cur_block_width, cur_block_height);
          save_buf_to_yuv(unfilter_block_buf, cur_block_height, cur_block_width,
                          unfilter_block_file);
          char filter_block_file[128];
          sprintf(filter_block_file,
                  "D:\\AOMedia\\aom_org\\test_result\\E19072501\\block_"
                  "yuv\\filter_block_%"
                  "d_%d_%dx%d.yuv",
                  x, y, cur_block_width, cur_block_height);
          save_buf_to_yuv(filter_block_buf, cur_block_height, cur_block_width,
                          filter_block_file);*/

          for (int i = 0; i < cur_block_height; i++) {
            for (int j = 0; j < cur_block_width; j++) {
              y_filter_buf[y * buf_height + i][x * buf_width + j] =
                  filter_block_buf[i][j];
#if out_yuv
              y_recon_buf[y * buf_height + i][x * buf_width + j] =
                  unfilter_block_buf[i][j];

              u_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                  *(recon_buf->u_buffer +
                    (i + y * buf_height) / 2 * recon_buf->uv_stride +
                    (j + x * buf_width) / 2);
              v_recon_buf[(y * buf_height + i) / 2][(x * buf_width + j) / 2] =
                  *(recon_buf->v_buffer +
                    (i + y * buf_height) / 2 * recon_buf->uv_stride +
                    (j + x * buf_width) / 2);
#endif
            }
          }
        }
      }
    }
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        *(bkuPy + j) = y_filter_buf[i][j];  // fill in the luma buffer again
      }
      bkuPy += dgr_stride;
    }
#if out_yuv
    char recon_file[256];
    sprintf(recon_file,
            "D:\\konglingyi\\work\\av1\\AV1\\example\\guangyao_20181127\\test_"
            "result\\E19081201"
            "\\FourPeople_1280x720_64x64_rdo.yuv");
    save_buf_to_yuv(y_recon_buf, height, width, recon_file,
                    cm->current_frame.frame_number);
    save_buf_to_yuv(u_recon_buf, height / 2, width / 2, recon_file, 1);
    save_buf_to_yuv(v_recon_buf, height / 2, width / 2, recon_file, 1);
#endif
  }

  // save_frame_yuv(filter_file, recon_buf, cm->current_frame.order_hint);
}

void addition_handle_frame_CRLC(YV12_BUFFER_CONFIG *source_frame,
                                AV1_COMMON *cm) {
  YV12_BUFFER_CONFIG *pcPicYuvRec = &cm->cur_frame->buf;

  if (!cm->seq_params.use_highbitdepth) {
    uint8_t *dgr = pcPicYuvRec->y_buffer;
    uint8_t *src = source_frame->y_buffer;

    int height = pcPicYuvRec->y_height;
    int width = pcPicYuvRec->y_width;
    int dgr_stride = pcPicYuvRec->y_stride;
    uint8_t **buf = new uint8_t *[height];
    for (int i = 0; i < height; i++) {
      buf[i] = new uint8_t[width];
    }
	
    buf = TF_Predict_CRLC(dgr, src, &cm->crlc_info[0], height, width, dgr_stride,
                          source_frame->y_stride, cm->base_qindex,
                          cm->crlc_info->crlc_unit_size);

    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        *(dgr + j) = buf[i][j];  // Fill in the luma buffer again
      }
      dgr += dgr_stride;
    }
  } else {
    uint16_t *py = CONVERT_TO_SHORTPTR(pcPicYuvRec->y_buffer);
    uint16_t *bkuPy = py;

    int height = pcPicYuvRec->y_height;
    int width = pcPicYuvRec->y_width;
    int stride = pcPicYuvRec->y_stride;

    uint16_t **buf = new uint16_t *[height];
    for (int i = 0; i < height; i++) {
      buf[i] = new uint16_t[width];
    }

    buf = TF_Predict_hbd(py, height, width, stride);

    // FILE *ff = fopen("end.yuv", "wb");
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        *(bkuPy + j) = buf[i][j];  // Fill in the luma buffer again
        // fwrite(bkuPy + j, 1, 1, ff);
      }
      bkuPy += stride;
    }
  }
  finish_python();
}
#endif

void block_rdo(YV12_BUFFER_CONFIG *src, YV12_BUFFER_CONFIG *filter,
               YV12_BUFFER_CONFIG *cnn) {
  int block_size = 64;
  int height = src->y_height;
  int width = src->y_width;
  for (int i = 0; i < height; i += block_size)
    for (int j = 0; j < width; j += block_size) {
      int rel_bs_h = i + block_size <= height ? block_size : height - i;
      int rel_bs_w = j + block_size <= width ? block_size : width - j;

      int64_t cnn_sse = aom_get_y_sse_part(src, cnn, i, rel_bs_h, j, rel_bs_w);
      int64_t filter_sse =
          aom_get_y_sse_part(src, filter, i, rel_bs_h, j, rel_bs_w);
      if (cnn_sse < filter_sse) {
        aom_yv12_partial_copy_y(cnn, i, i + rel_bs_h, j, j + rel_bs_w, filter,
                                i, j);
      } else {
	  }
    }
}