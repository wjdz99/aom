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
set(AOM_ARCH "mips32")
set(ARCH_MIPS 1)

# Assembly flavors available for the target.
set(HAVE_MIPS32 1)
set(HAVE_DSPR2 1)

# RTCD versions of assembly flavor flags ("yes" means on in rtcd.pl, not 1).
set(RTCD_ARCH_MIPS "yes")
set(RTCD_HAVE_MIPS32 "yes")
set(RTCD_HAVE_DSPR2 "yes")

foreach (config_var ${AOM_CONFIG_VARS})
  if (${${config_var}})
    set(RTCD_${config_var} yes)
  endif ()
endforeach ()
