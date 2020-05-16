#!/usr/bin/env python
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
__author__ = "maggie.sun@intel.com, ryan.lei@intel.com"

import os
import Utils
from Config import BinPath, LogCmdOnly
from Utils import ExecuteCmd

#use ffmpeg to decode bitstream
def DecodeWithFfmpeg(infile, outfile):
    decoder = os.path.join(BinPath, 'ffmpeg.exe')
    args = " -y -i %s -pix_fmt yuv420p -c:v rawvideo %s" % (infile, outfile)
    cmd = decoder + args
    ExecuteCmd(cmd, LogCmdOnly)

def DecodeWithAOM(infile, outfile):
    decoder = os.path.join(BinPath, 'aomdec.exe')
    args = " --codec=av1 --i420 --rawvideo --summary -o %s %s" % (outfile, infile)
    cmd = decoder + args
    ExecuteCmd(cmd, LogCmdOnly)

def VideoDecode(codec, infile, outfile):
    Utils.CmdLogger.write("::Decode\n")
    if codec == 'AV1':
        DecodeWithAOM(infile, outfile)
    else:
        DecodeWithFfmpeg(infile, outfile)

