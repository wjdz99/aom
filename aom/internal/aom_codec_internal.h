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

/*!\file
 * \brief Describes the decoder algorithm interface for algorithm
 *        implementations.
 *
 * This file defines the private structures and data types that are only
 * relevant to implementing an algorithm, as opposed to using it.
 *
 * To create a decoder algorithm class, an interface structure is put
 * into the global namespace:
 *     <pre>
 *     my_codec.c:
 *       AomCodecIfaceT my_codec = {
 *           "My Codec v1.0",
 *           AOM_CODEC_ALG_ABI_VERSION,
 *           ...
 *       };
 *     </pre>
 *
 * An application instantiates a specific decoder instance by using
 * aom_codec_init() and a pointer to the algorithm's interface structure:
 *     <pre>
 *     my_app.c:
 *       extern AomCodecIfaceT my_codec;
 *       {
 *           AomCodecCtxT algo;
 *           res = aom_codec_init(&algo, &my_codec);
 *       }
 *     </pre>
 *
 * Once initialized, the instance is manged using other functions from
 * the aom_codec_* family.
 */
#ifndef AOM_INTERNAL_AOM_CODEC_INTERNAL_H_
#define AOM_INTERNAL_AOM_CODEC_INTERNAL_H_
#include "./aom_config.h"
#include "../aom_decoder.h"
#include "../aom_encoder.h"
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!\brief Current ABI version number
 *
 * \internal
 * If this file is altered in any way that changes the ABI, this value
 * must be bumped.  Examples include, but are not limited to, changing
 * types, removing or reassigning enums, adding/removing/rearranging
 * fields to structures
 */
#define AOM_CODEC_INTERNAL_ABI_VERSION (5) /**<\hideinitializer*/

typedef struct aom_codec_alg_priv AomCodecAlgPrivT;
typedef struct AomCodecPrivEncMrCfg AomCodecPrivEncMrCfgT;

/*!\brief init function pointer prototype
 *
 * Performs algorithm-specific initialization of the decoder context. This
 * function is called by the generic aom_codec_init() wrapper function, so
 * plugins implementing this interface may trust the input parameters to be
 * properly initialized.
 *
 * \param[in] ctx   Pointer to this instance's context
 * \retval #AOM_CODEC_OK
 *     The input stream was recognized and decoder initialized.
 * \retval #AOM_CODEC_MEM_ERROR
 *     Memory operation failed.
 */
typedef AomCodecErrT (*AomCodecInitFnT)(AomCodecCtxT *ctx,
                                        AomCodecPrivEncMrCfgT *data);

/*!\brief destroy function pointer prototype
 *
 * Performs algorithm-specific destruction of the decoder context. This
 * function is called by the generic aom_codec_destroy() wrapper function,
 * so plugins implementing this interface may trust the input parameters
 * to be properly initialized.
 *
 * \param[in] ctx   Pointer to this instance's context
 * \retval #AOM_CODEC_OK
 *     The input stream was recognized and decoder initialized.
 * \retval #AOM_CODEC_MEM_ERROR
 *     Memory operation failed.
 */
typedef AomCodecErrT (*AomCodecDestroyFnT)(AomCodecAlgPrivT *ctx);

/*!\brief parse stream info function pointer prototype
 *
 * Performs high level parsing of the bitstream. This function is called by the
 * generic aom_codec_peek_stream_info() wrapper function, so plugins
 * implementing this interface may trust the input parameters to be properly
 * initialized.
 *
 * \param[in]      data    Pointer to a block of data to parse
 * \param[in]      data_sz Size of the data buffer
 * \param[in,out]  si      Pointer to stream info to update. The size member
 *                         \ref MUST be properly initialized, but \ref MAY be
 *                         clobbered by the algorithm. This parameter \ref MAY
 *                         be NULL.
 *
 * \retval #AOM_CODEC_OK
 *     Bitstream is parsable and stream information updated
 */
typedef AomCodecErrT (*AomCodecPeekSiFnT)(const uint8_t *data,
                                          unsigned int data_sz,
                                          AomCodecStreamInfoT *si);

/*!\brief Return information about the current stream.
 *
 * Returns information about the stream that has been parsed during decoding.
 *
 * \param[in]      ctx     Pointer to this instance's context
 * \param[in,out]  si      Pointer to stream info to update. The size member
 *                         \ref MUST be properly initialized, but \ref MAY be
 *                         clobbered by the algorithm. This parameter \ref MAY
 *                         be NULL.
 *
 * \retval #AOM_CODEC_OK
 *     Bitstream is parsable and stream information updated
 */
typedef AomCodecErrT (*AomCodecGetSiFnT)(AomCodecAlgPrivT *ctx,
                                         AomCodecStreamInfoT *si);

/*!\brief control function pointer prototype
 *
 * This function is used to exchange algorithm specific data with the decoder
 * instance. This can be used to implement features specific to a particular
 * algorithm.
 *
 * This function is called by the generic aom_codec_control() wrapper
 * function, so plugins implementing this interface may trust the input
 * parameters to be properly initialized. However,  this interface does not
 * provide type safety for the exchanged data or assign meanings to the
 * control codes. Those details should be specified in the algorithm's
 * header file. In particular, the ctrl_id parameter is guaranteed to exist
 * in the algorithm's control Mapping table, and the data parameter may be NULL.
 *
 *
 * \param[in]     ctx              Pointer to this instance's context
 * \param[in]     ctrl_id          Algorithm specific control identifier
 * \param[in,out] data             Data to exchange with algorithm instance.
 *
 * \retval #AOM_CODEC_OK
 *     The internal state data was deserialized.
 */
typedef AomCodecErrT (*AomCodecControlFnT)(AomCodecAlgPrivT *ctx, va_list ap);

/*!\brief control function pointer Mapping
 *
 * This structure stores the Mapping between control identifiers and
 * implementing functions. Each algorithm provides a list of these
 * mappings. This list is searched by the aom_codec_control() wrapper
 * function to determine which function to invoke. The special
 * value {0, NULL} is used to indicate end-of-list, and must be
 * present. The special value {0, <non-null>} can be used as a catch-all
 * Mapping. This implies that ctrl_id values chosen by the algorithm
 * \ref MUST be non-zero.
 */
typedef const struct AomCodecCtrlFnMap {
  int ctrl_id;
  AomCodecControlFnT fn;
} AomCodecCtrlFnMapT;

/*!\brief decode data function pointer prototype
 *
 * Processes a buffer of coded data. If the processing results in a new
 * decoded frame becoming available, #AOM_CODEC_CB_PUT_SLICE and
 * #AOM_CODEC_CB_PUT_FRAME events are generated as appropriate. This
 * function is called by the generic aom_codec_decode() wrapper function,
 * so plugins implementing this interface may trust the input parameters
 * to be properly initialized.
 *
 * \param[in] ctx          Pointer to this instance's context
 * \param[in] data         Pointer to this block of new coded data. If
 *                         NULL, a #AOM_CODEC_CB_PUT_FRAME event is posted
 *                         for the previously decoded frame.
 * \param[in] data_sz      Size of the coded data, in bytes.
 *
 * \return Returns #AOM_CODEC_OK if the coded data was processed completely
 *         and future pictures can be decoded without error. Otherwise,
 *         see the descriptions of the other error codes in ::AomCodecErrT
 *         for recoverability capabilities.
 */
typedef AomCodecErrT (*AomCodecDecodeFnT)(AomCodecAlgPrivT *ctx,
                                          const uint8_t *data,
                                          unsigned int data_sz, void *user_priv,
                                          long deadline);

/*!\brief Decoded frames iterator
 *
 * Iterates over a list of the frames available for display. The iterator
 * storage should be initialized to NULL to start the iteration. Iteration is
 * complete when this function returns NULL.
 *
 * The list of available frames becomes valid upon completion of the
 * aom_codec_decode call, and remains valid until the next call to
 * aom_codec_decode.
 *
 * \param[in]     ctx      Pointer to this instance's context
 * \param[in out] iter     Iterator storage, initialized to NULL
 *
 * \return Returns a pointer to an image, if one is ready for display. Frames
 *         produced will always be in PTS (presentation time stamp) order.
 */
typedef AomImageT *(*AomCodecGetFrameFnT)(AomCodecAlgPrivT *ctx,
                                          AomCodecIterT *iter);

/*!\brief Pass in external frame buffers for the decoder to use.
 *
 * Registers functions to be called when libaom needs a frame buffer
 * to decode the current frame and a function to be called when libaom does
 * not internally reference the frame buffer. This set function must
 * be called before the first call to decode or libaom will assume the
 * default behavior of allocating frame buffers internally.
 *
 * \param[in] ctx          Pointer to this instance's context
 * \param[in] cb_get       Pointer to the get callback function
 * \param[in] cb_release   Pointer to the release callback function
 * \param[in] cb_priv      Callback's private data
 *
 * \retval #AOM_CODEC_OK
 *     External frame buffers will be used by libaom.
 * \retval #AOM_CODEC_INVALID_PARAM
 *     One or more of the callbacks were NULL.
 * \retval #AOM_CODEC_ERROR
 *     Decoder context not initialized, or algorithm not capable of
 *     using external frame buffers.
 *
 * \note
 * When decoding AV1, the application may be required to pass in at least
 * #AOM_MAXIMUM_WORK_BUFFERS external frame
 * buffers.
 */
typedef AomCodecErrT (*AomCodecSetFbFnT)(AomCodecAlgPrivT *ctx,
                                         AomGetFrameBufferCbFnT cb_get,
                                         AomReleaseFrameBufferCbFnT cb_release,
                                         void *cb_priv);

typedef AomCodecErrT (*AomCodecEncodeFnT)(
    AomCodecAlgPrivT *ctx, const AomImageT *img, AomCodecPtsT pts,
    unsigned long duration, AomEncFrameFlagsT flags, unsigned long deadline);
typedef const AomCodecCxPktT *(*AomCodecGetCxDataFnT)(AomCodecAlgPrivT *ctx,
                                                      AomCodecIterT *iter);

typedef AomCodecErrT (*AomCodecEncConfigSetFnT)(AomCodecAlgPrivT *ctx,
                                                const AomCodecEncCfgT *cfg);
typedef AomFixedBufT *(*AomCodecGetGlobalHeadersFnT)(AomCodecAlgPrivT *ctx);

typedef AomImageT *(*AomCodecGetPreviewFrameFnT)(AomCodecAlgPrivT *ctx);

typedef AomCodecErrT (*AomCodecEncMrGetMemLocFnT)(const AomCodecEncCfgT *cfg,
                                                  void **mem_loc);

/*!\brief usage configuration Mapping
 *
 * This structure stores the Mapping between usage identifiers and
 * configuration structures. Each algorithm provides a list of these
 * mappings. This list is searched by the aom_codec_enc_config_default()
 * wrapper function to determine which config to return. The special value
 * {-1, {0}} is used to indicate end-of-list, and must be present. At least
 * one Mapping must be present, in addition to the end-of-list.
 *
 */
typedef const struct AomCodecEncCfgMap {
  int usage;
  AomCodecEncCfgT cfg;
} AomCodecEncCfgMapT;

/*!\brief Decoder algorithm interface interface
 *
 * All decoders \ref MUST expose a variable of this type.
 */
struct AomCodecIface {
  const char *name;              /**< Identification String  */
  int abi_version;               /**< Implemented ABI version */
  AomCodecCapsT caps;            /**< Decoder capabilities */
  AomCodecInitFnT init;          /**< \copydoc ::AomCodecInitFnT */
  AomCodecDestroyFnT destroy;    /**< \copydoc ::AomCodecDestroyFnT */
  AomCodecCtrlFnMapT *ctrl_maps; /**< \copydoc ::AomCodecCtrlFnMapT */
  struct AomCodecDecIface {
    AomCodecPeekSiFnT peek_si;     /**< \copydoc ::AomCodecPeekSiFnT */
    AomCodecGetSiFnT get_si;       /**< \copydoc ::AomCodecGetSiFnT */
    AomCodecDecodeFnT decode;      /**< \copydoc ::AomCodecDecodeFnT */
    AomCodecGetFrameFnT get_frame; /**< \copydoc ::AomCodecGetFrameFnT */
    AomCodecSetFbFnT set_fb_fn;    /**< \copydoc ::AomCodecSetFbFnT */
  } dec;
  struct AomCodecEncIface {
    int cfg_map_count;
    AomCodecEncCfgMapT *cfg_maps;     /**< \copydoc ::AomCodecEncCfgMapT */
    AomCodecEncodeFnT encode;         /**< \copydoc ::AomCodecEncodeFnT */
    AomCodecGetCxDataFnT get_cx_data; /**< \copydoc ::AomCodecGetCxDataFnT */
    AomCodecEncConfigSetFnT cfg_set;  /**< \copydoc ::AomCodecEncConfigSetFnT */
    AomCodecGetGlobalHeadersFnT
        get_glob_hdrs; /**< \copydoc ::AomCodecGetGlobalHeadersFnT */
    AomCodecGetPreviewFrameFnT
        get_preview; /**< \copydoc ::AomCodecGetPreviewFrameFnT */
    AomCodecEncMrGetMemLocFnT
        mr_get_mem_loc; /**< \copydoc ::AomCodecEncMrGetMemLocFnT */
  } enc;
};

/*!\brief Callback function pointer / user data pair storage */
typedef struct AomCodecPrivCbPair {
  union {
    AomCodecPutFrameCbFnT put_frame;
    AomCodecPutSliceCbFnT put_slice;
  } u;
  void *user_priv;
} AomCodecPrivCbPairT;

/*!\brief Instance private storage
 *
 * This structure is allocated by the algorithm's init function. It can be
 * extended in one of two ways. First, a second, algorithm specific structure
 * can be allocated and the priv member pointed to it. Alternatively, this
 * structure can be made the first member of the algorithm specific structure,
 * and the pointer cast to the proper type.
 */
struct AomCodecPriv {
  const char *err_detail;
  AomCodecFlagsT init_flags;
  struct {
    AomCodecPrivCbPairT put_frame_cb;
    AomCodecPrivCbPairT put_slice_cb;
  } dec;
  struct {
    AomFixedBufT cx_data_dst_buf;
    unsigned int cx_data_pad_before;
    unsigned int cx_data_pad_after;
    AomCodecCxPktT cx_data_pkt;
    unsigned int total_encoders;
  } enc;
};

/*
 * Multi-resolution encoding internal configuration
 */
struct AomCodecPrivEncMrCfg {
  unsigned int mr_total_resolutions;
  unsigned int mr_encoder_id;
  struct AomRational mr_down_sampling_factor;
  void *mr_low_res_mode_info;
};

#undef AOM_CTRL_USE_TYPE
#define AOM_CTRL_USE_TYPE(id, typ) \
  static AOM_INLINE typ id##__value(va_list args) { return va_arg(args, typ); }

#undef AOM_CTRL_USE_TYPE_DEPRECATED
#define AOM_CTRL_USE_TYPE_DEPRECATED(id, typ) \
  static AOM_INLINE typ id##__value(va_list args) { return va_arg(args, typ); }

#define CAST(id, arg) id##__value(arg)

/* CODEC_INTERFACE convenience macro
 *
 * By convention, each codec interface is a struct with extern linkage, where
 * the symbol is suffixed with _algo. A getter function is also defined to
 * return a pointer to the struct, since in some cases it's easier to work
 * with text symbols than data symbols (see issue #169). This function has
 * the same name as the struct, less the _algo suffix. The CODEC_INTERFACE
 * macro is provided to define this getter function automatically.
 */
#define CODEC_INTERFACE(id)                       \
  AomCodecIfaceT *id(void) { return &id##_algo; } \
  AomCodecIfaceT id##_algo

/* Internal Utility Functions
 *
 * The following functions are intended to be used inside algorithms as
 * utilities for manipulating aom_codec_* data structures.
 */
struct AomCodecPktList {
  unsigned int cnt;
  unsigned int max;
  struct AomCodecCxPkt pkts[1];
};

#define aom_codec_pkt_list_decl(n)  \
  union {                           \
    struct AomCodecPktList head;    \
    struct {                        \
      struct AomCodecPktList head;  \
      struct AomCodecCxPkt pkts[n]; \
    } alloc;                        \
  }

#define aom_codec_pkt_list_init(m) \
  (m)->alloc.head.cnt = 0,         \
  (m)->alloc.head.max = sizeof((m)->alloc.pkts) / sizeof((m)->alloc.pkts[0])

int aom_codec_pkt_list_add(struct AomCodecPktList *,
                           const struct AomCodecCxPkt *);

const AomCodecCxPktT *aom_codec_pkt_list_get(struct AomCodecPktList *list,
                                             AomCodecIterT *iter);

#include <stdio.h>
#include <setjmp.h>

struct AomInternalErrorInfo {
  AomCodecErrT error_code;
  int has_detail;
  char detail[80];
  int setjmp;
  jmp_buf jmp;
};

#define CLANG_ANALYZER_NORETURN
#if defined(__has_feature)
#if __has_feature(attribute_analyzer_noreturn)
#undef CLANG_ANALYZER_NORETURN
#define CLANG_ANALYZER_NORETURN __attribute__((analyzer_noreturn))
#endif
#endif

void aom_internal_error(struct AomInternalErrorInfo *info, AomCodecErrT error,
                        const char *fmt, ...) CLANG_ANALYZER_NORETURN;

void aom_merge_corrupted_flag(int *corrupted, int value);

#if CONFIG_DEBUG
#define AOM_CHECK_MEM_ERROR(error_info, lval, expr)                         \
  do {                                                                      \
    lval = (expr);                                                          \
    if (!lval)                                                              \
      aom_internal_error(error_info, AOM_CODEC_MEM_ERROR,                   \
                         "Failed to allocate " #lval " at %s:%d", __FILE__, \
                         __LINE__);                                         \
  } while (0)
#else
#define AOM_CHECK_MEM_ERROR(error_info, lval, expr)       \
  do {                                                    \
    lval = (expr);                                        \
    if (!lval)                                            \
      aom_internal_error(error_info, AOM_CODEC_MEM_ERROR, \
                         "Failed to allocate " #lval);    \
  } while (0)
#endif
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_INTERNAL_AOM_CODEC_INTERNAL_H_
