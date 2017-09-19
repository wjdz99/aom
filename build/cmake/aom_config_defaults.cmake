##
## Copyright (c) 2016, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
# Defaults for every libaom configuration variable. Here we add all libaom
# config variables to the cmake variable cache, but omit the FORCE parameter to
# allow users to specify values when executing cmake to generate build files.
# Values here are used only if not set by the user.
set(INLINE "" CACHE STRING "Sets INLINE value for current target.")
set(RESTRICT "" CACHE STRING "Sets RESTRICT value for current target.")

# CPUs.
set(ARCH_ARM 0 CACHE NUMBER "Enables ARM architecture.")
set(ARCH_MIPS 0 CACHE NUMBER "Enables MIPS architecture.")
set(ARCH_X86 0 CACHE NUMBER "Enables X86 architecture.")
set(ARCH_X86_64 0 CACHE NUMBER "Enables X86_64 architecture.")

# ARM optimization flags.
set(HAVE_NEON 0 CACHE NUMBER "Enables NEON intrinsics optimizations.")
set(HAVE_NEON_ASM 0 CACHE NUMBER "Enables NEON assembly optimizations.")

# MIPS optimization flags.
set(HAVE_DSPR2 0 CACHE NUMBER "Enables DSPR2 optimizations.")
set(HAVE_MIPS32 0 CACHE NUMBER "Enables MIPS32 optimizations.")
set(HAVE_MIPS64 0 CACHE NUMBER "Enables MIPS64 optimizations. ")
set(HAVE_MSA 0 CACHE NUMBER "Enables MSA optimizations.")

# x86/x86_64 optimization flags.
set(HAVE_AVX 0 CACHE NUMBER "Enables AVX optimizations.")
set(HAVE_AVX2 0 CACHE NUMBER "Enables AVX2 optimizations.")
set(HAVE_MMX 0 CACHE NUMBER "Enables MMX optimizations. ")
set(HAVE_SSE 0 CACHE NUMBER "Enables SSE optimizations.")
set(HAVE_SSE2 0 CACHE NUMBER "Enables SSE2 optimizations.")
set(HAVE_SSE3 0 CACHE NUMBER "Enables SSE3 optimizations.")
set(HAVE_SSE4_1 0 CACHE NUMBER "Enables SSE 4.1 optimizations.")
set(HAVE_SSSE3 0 CACHE NUMBER "Enables SSSE3 optimizations.")

# Flags describing the build environment.
set(HAVE_AOM_PORTS 0 CACHE NUMBER "Internal flag, deprecated.")
set(HAVE_FEXCEPT 0 CACHE NUMBER "Internal flag, GNU fenv.h present for target.")
set(HAVE_PTHREAD_H 0 CACHE NUMBER "Internal flag, target pthread support.")
set(HAVE_UNISTD_H 0 CACHE NUMBER "Internal flag, unistd.h present for target.")
set(HAVE_WXWIDGETS 0 CACHE NUMBER "WxWidgets present.")

# Deprecated flags preserved for compatibility with configure build.
set(CONFIG_CODEC_SRCS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_DEBUG_LIBS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_DEPENDENCY_TRACKING 1 CACHE NUMBER "Deprecated flag.")
set(CONFIG_EXPERIMENTAL 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_EXTERNAL_BUILD 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_INSTALL_BINS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_INSTALL_DOCS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_INSTALL_LIBS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_INSTALL_SRCS 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_POSTPROC 1 CACHE NUMBER "Deprecated flag.")
set(CONFIG_POSTPROC_VISUALIZER 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_RVCT 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_SMALL 0 CACHE NUMBER "Deprecated flag.")
set(CONFIG_STATIC_MSVCRT 0 CACHE NUMBER "Deprecated flag.")

# Build configuration flags.
set(CONFIG_AV1 1 CACHE NUMBER "Internal flag.")
set(CONFIG_AV1_DECODER 1 CACHE NUMBER "Enable AV1 decoder.")
set(CONFIG_AV1_ENCODER 1 CACHE NUMBER "Enable AV1 encoder.")
set(CONFIG_BIG_ENDIAN 0 CACHE NUMBER "Internal flag.")
set(CONFIG_GCC 0 CACHE NUMBER "Building with GCC (detected).")
set(CONFIG_GCOV 0 CACHE NUMBER "Enable gcov support.")
set(CONFIG_GPROF 0 CACHE NUMBER "Enable gprof support.")
set(CONFIG_LIBYUV 1 CACHE NUMBER "Enables libyuv scaling/conversion support.")
set(CONFIG_MSVS 0 CACHE NUMBER "Building with MS Visual Studio (detected).")
set(CONFIG_MULTITHREAD 1 CACHE NUMBER "Multithread support.")
set(CONFIG_OS_SUPPORT 0 CACHE NUMBER "Internal flag.")
set(CONFIG_PIC 0 CACHE NUMBER "Build with PIC enabled.")
set(CONFIG_RUNTIME_CPU_DETECT 1 CACHE NUMBER "Runtime CPU detection support.")
set(CONFIG_SHARED 0 CACHE NUMBER "Build shared libs.")
set(CONFIG_STATIC 1 CACHE NUMBER "Build static libs.")
set(CONFIG_WEBM_IO 1 CACHE NUMBER "Enables WebM support.")

# Debugging flags.
set(CONFIG_BITSTREAM_DEBUG 0 CACHE NUMBER "Bitstream debugging flag.")
set(CONFIG_DEBUG 0 CACHE NUMBER "Debug build flag.")

# Testing flags.
set(CONFIG_DECODE_PERF_TESTS 0 CACHE NUMBER "Enables decoder performance test.")
set(CONFIG_ENCODE_PERF_TESTS 0 CACHE NUMBER "Enables encoder performance test.")
set(CONFIG_UNIT_TESTS 1 CACHE NUMBER "Enables unit tests.")

# AV1 feature flags.
set(CONFIG_ACCOUNTING 0 CACHE NUMBER "Enables bit accounting.")
set(CONFIG_ANALYZER 0 CACHE NUMBER "Enables bit stream analyzer.")
set(CONFIG_COEFFICIENT_RANGE_CHECKING 0 CACHE NUMBER "Coefficient range check.")
set(CONFIG_HIGHBITDEPTH 1 CACHE NUMBER "Enables high bit depth support.")
set(CONFIG_INSPECTION 0 CACHE NUMBER "Enables bitstream inspection.")
set(CONFIG_INTERNAL_STATS 0 CACHE NUMBER "Codec stats.")
set(CONFIG_LOWBITDEPTH 1 CACHE NUMBER "Enables low bit depth support.")
set(CONFIG_REALTIME_ONLY 0 CACHE NUMBER "Support only realtime encodes.")
set(CONFIG_SIZE_LIMIT 0 CACHE NUMBER "Limit max decode width/height.")
set(CONFIG_SPATIAL_RESAMPLING 1 CACHE NUMBER "Spatial resampling.")

# AV1 experiment flags.
set(CONFIG_ADAPT_SCAN 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_AMVR 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_ANS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_AOM_QM 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_BGSPRITE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CB4X4 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CDEF 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CDEF_SINGLEPASS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CFL 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CHROMA_2X2 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CHROMA_SUB8X8 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_COEF_INTERLEAVE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_COLORSPACE_HEADERS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_COMPOUND_ROUND 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_COMPOUND_SEGMENT 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_COMPOUND_SINGLEREF 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CONVOLVE_ROUND 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_CTX1D 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_DCT16 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_DCT32 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_DCT4 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_DCT64 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_DCT8 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DAALA_TX 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DCT_ONLY 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DEBLOCK_13TAP 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DEPENDENT_HORZTILES 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DIST_8X8 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DPCM_INTRA 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_DUAL_FILTER 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_ENTROPY_STATS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_COMP_REFS 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_DELTA_Q 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_INTER 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_INTRA 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_PARTITION 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_PARTITION_TYPES 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_PARTITION_TYPES_AB 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_REFS 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_TILE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_EXT_TX 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_FILTER_INTRA 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_FP_MB_STATS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_FRAME_MARKER 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_FRAME_SIZE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_FRAME_SUPERRES 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_GF_GROUPS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_GLOBAL_MOTION 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_HASH_ME 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_INTERINTRA 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_INTER_STATS_ONLY 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_INTRABC 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_INTRA_EDGE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_INTRA_INTERP 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_JNT_COMP 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LGT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LOOPFILTERING_ACROSS_TILES 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LOOPFILTER_LEVEL 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LOOP_RESTORATION 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LPF_DIRECT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LPF_SB 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_LV_MAP 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MASKED_TX 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MAX_TILE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MFMV 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MOTION_VAR 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MRC_TX 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_MV_COMPRESS 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_NCOBMC 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_NCOBMC_ADAPT_WEIGHT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_NEW_MULTISYMBOL 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_NEW_QUANT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_NO_FRAME_CONTEXT_SIGNALING 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_OBU 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_ONE_SIDED_COMPOUND 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_OPT_REF_MV 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_PALETTE_DELTA_ENCODING 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_PALETTE_THROUGHPUT 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_PARALLEL_DEBLOCKING 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_PVQ 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_Q_ADAPT_PROBS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_RAWBITS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_RD_DEBUG 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_RECT_TX 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_RECT_TX_EXT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_REFERENCE_BUFFER 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_REF_ADAPT 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_RESTRICT_COMPRESSED_HDR 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_SBL_SYMBOL 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_SEGMENT_ZEROMV 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_SMOOTH_HV 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_STRIPED_LOOP_RESTORATION 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_SUPERTX 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TEMPMV_SIGNALING 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TPL_MV 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TRIPRED 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TX64X64 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TXK_SEL 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_TXMG 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_UNPOISON_PARTITION_CTX 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_VAR_REFS 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_VAR_TX 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_VAR_TX_NO_TX_MODE 0 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_WARPED_MOTION 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_WEDGE 1 CACHE NUMBER "AV1 experiment flag.")
set(CONFIG_XIPHRC 0 CACHE NUMBER "AV1 experiment flag.")
