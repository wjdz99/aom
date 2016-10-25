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

cmake_minimum_required(VERSION 3.2)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

function (add_flag_if_supported flag)
  unset(C_FLAG_SUPPORTED CACHE)
  unset(CXX_FLAG_SUPPORTED CACHE)
  message("Checking C compiler flag support for: " ${flag})
  CHECK_C_COMPILER_FLAG("${flag}" C_FLAG_SUPPORTED)
  message("Checking CXX compiler flag support for: " ${flag})
  CHECK_CXX_COMPILER_FLAG("${flag}" CXX_FLAG_SUPPORTED)
  if (C_FLAG_SUPPORTED AND CXX_FLAG_SUPPORTED)
    set(CMAKE_C_FLAGS "${flag} ${CMAKE_C_FLAGS}" CACHE STRING "" FORCE)
    set(CMAKE_CXX_FLAGS "${flag} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
  endif ()
endfunction ()

# Set warning levels.
if (MSVC)
  add_flag_if_supported("/W4")
  # Disable MSVC warnings that suggest making code non-portable.
  add_flag_if_supported("/wd4996")
  if (ENABLE_WERROR)
    add_flag_if_supported("/WX")
  endif ()
else ()
  add_flag_if_supported("-Wall")
  add_flag_if_supported("-Wextra")
  add_flag_if_supported("-Wno-deprecated")
  add_flag_if_supported("-Wshorten-64-to-32")
  add_flag_if_supported("-Wnarrowing")
  if (ENABLE_WERROR)
    add_flag_if_supported("-Werror")
  endif ()
endif ()

