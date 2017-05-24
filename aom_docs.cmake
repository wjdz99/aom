##
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
if (NOT AOM_AOM_DOCS_CMAKE_)
set(AOM_AOM_DOCS_CMAKE_ 1)

cmake_minimum_required(VERSION 3.5)

include(FindDoxygen)
if (NOT DOXYGEN_FOUND)
  message(FATAL_ERROR "Can't find doxygen, install it or disable ENABLE_DOCS")
endif ()

set(AOM_DOXYGEN_CONFIG_TEMPLATE "${AOM_ROOT}/libs.doxy_template")
set(AOM_DOXYGEN_SOURCES
    "${AOM_ROOT}/aom/aom.h"
    "${AOM_ROOT}/keywords.dox"
    "${AOM_ROOT}/mainpage.dox"
    "${AOM_ROOT}/usage.dox")

if (CONFIG_AV1_DECODER)
  set(AOM_DOXYGEN_EXAMPLE_SOURCES
      ${AOM_DOXYGEN_EXAMPLE_SOURCES}
      "${AOM_ROOT}/aomdec.c"
      "${AOM_ROOT}/examples/decode_to_md5.c"
      "${AOM_ROOT}/examples/decode_with_drops.c"
      "${AOM_ROOT}/examples/simple_decoder.c")

  set(AOM_DOXYGEN_SOURCES
      ${AOM_DOXYGEN_SOURCES}
      "${AOM_ROOT}/aom/aomdx.h"
      "${AOM_ROOT}/usage_dx.dox")

  if (CONFIG_ANALYZER)
    set(AOM_DOXYGEN_EXAMPLE_SOURCES 
        ${AOM_DOXYGEN_EXAMPLE_SOURCES}
        "${AOM_ROOT}examples/analyzer.cc")
  endif ()

  if (CONFIG_INSPECTION)
     set(AOM_DOXYGEN_EXAMPLE_SOURCES
         ${AOM_DOXYGEN_EXAMPLE_SOURCES}
         "${AOM_ROOT}/examples/inspect.c")
  endif ()
endif ()

if (CONFIG_AV1_ENCODER)
  set(AOM_DOXYGEN_SOURCES
      ${AOM_DOXYGEN_SOURCES}
      "${AOM_ROOT}/aom/aomcx.h"
      "${AOM_ROOT}/usage_cx.dox")

  set(AOM_DOXYGEN_EXAMPLE_SOURCES
      ${AOM_DOXYGEN_EXAMPLE_SOURCES}
      "${AOM_ROOT}/aomenc.c"
      "${AOM_ROOT}/examples/lossless_encoder.c"
      "${AOM_ROOT}/examples/set_maps.c"
      "${AOM_ROOT}/examples/simple_encoder.c"
      "${AOM_ROOT}/examples/twopass_encoder.c")
endif ()


endif ()

if (CONFIG_AV1_DECODER AND CONFIG_AV1_ENCODER)
  set(AOM_DOXYGEN_EXAMPLE_SOURCES
      ${AOM_DOXYGEN_EXAMPLE_SOURCES}
      "${AOM_ROOT}/examples/aom_cx_set_ref.c")
endif ()

function (setup_doxygen_targets)
  if (CONFIG_AV1_DECODER)
  endif ()
  if (CONFIG_AV1_ENCODER)
  endif ()
  add_custom_target(generate_docs)
endfunction ()

endif ()  # AOM_AOM_DOCS_CMAKE_
