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

# Adds an object library target. Terminates generation if $flag is not supported
# by the current compiler. $flag is the intrinsics flag required by the current
# compiler, and is added to the compile flags for all sources in $sources.
# $opt_name is used to name the target. $target_to_update is made
# dependent upon the created target.
function (add_intrinsics_target flag opt_name target_to_update sources)
  require_flag(${flag})
  set(target_name aom_${opt_name}_intrinsics)
  add_library(${target_name} OBJECT ${${sources}})
  target_compile_options(${target_name} PUBLIC ${flag})
  target_sources(${target_to_update} PUBLIC $<TARGET_OBJECTS:${target_name}>)
endfunction ()
