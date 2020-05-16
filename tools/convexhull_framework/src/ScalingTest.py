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
import xlsxwriter
import xlrd
import logging
from CalculateQualityMetrics import CalculateQualityMetric, GatherQualityMetrics
from VideoScaler import GetDownScaledOutFile, DownScaling, UpScaling, GetUpScaledOutFile
from Config import DnScaleRatio, FrameNum, Clips, QualityList, LoggerName, Path_ScalingResults, ScalQty_WtCols, \
    ScalQty_startRow, LineColors, ScalSumQty_WtCols, ContentPath
from Utils import GetVideoInfo, Cleanfolder, GetShortContentName, CreateChart_Scatter, AddSeriesToChart_Scatter, UpdateChart, \
    InsertChartsToSheet, CalcRowsClassAndContentDict, GetContents, CreateChart_Line, AddSeriesToChart_Line

subloggername = "ScalingTest"
loggername = LoggerName + '.' + '%s' % subloggername
logger = logging.getLogger(loggername)

StatsMetrics = ['Max', 'Std+', 'Average', 'Std-', 'Min']

def GetScalingResultExcelFile(rationum, scalAlgoNum):
    filename = "ScalingResults_RatioNum_%d_AlgoNum_%d.xlsx" % (rationum, scalAlgoNum)
    file = os.path.join(Path_ScalingResults, filename)
    return file

def GetScalingResultExcelFile_PerContent(content):
    filename = GetShortContentName(content)
    filename = "ScalingResults_%s.xlsx" % filename
    file = os.path.join(Path_ScalingResults, filename)
    return file

def Run_Scaling_Test(content, dnScalAlgo, upScalAlgo, path_dnscl, path_upscl, path_log, path_cfg,
                     savememory, keepupscaledyuv):
    logger.info("start running scaling tests with content %s" % os.path.basename(content))
    cls, width, height, fr, bitdepth, fmt, totalnum = GetVideoInfo(content, Clips)
    if totalnum < FrameNum:
        logger.error("content %s includes total %d frames < specified % frames!" % (content, totalnum, FrameNum))
        return

    DnScaledRes = [(int(width / ratio), int(height / ratio)) for ratio in DnScaleRatio]
    for i in range(len(DnScaledRes)):
        if savememory:
            Cleanfolder(path_cfg)
            if not keepupscaledyuv:
                Cleanfolder(path_upscl)

        DnScaledW = DnScaledRes[i][0]
        DnScaledH = DnScaledRes[i][1]
        logger.info("start downscaling content to %dx%d" % (DnScaledW, DnScaledH))
        # downscaling
        dnscalyuv = GetDownScaledOutFile(content, width, height, DnScaledW, DnScaledH, path_dnscl, dnScalAlgo)
        if not os.path.isfile(dnscalyuv):
            dnscalyuv = DownScaling(content, FrameNum, width, height, DnScaledW, DnScaledH, path_dnscl, dnScalAlgo)

        upscaleyuv = UpScaling(dnscalyuv, FrameNum, DnScaledW, DnScaledH, width, height, path_upscl, upScalAlgo)

        CalculateQualityMetric(content, FrameNum, upscaleyuv, width, height, path_log, path_cfg)

    if savememory:
        Cleanfolder(path_dnscl)

    logger.info("finish running scaling test.")

def GeneratePerContentExcelFile(dnScalAlgos, upScalAlgos, content, scaleRatios, path_log):
    contshortname = GetShortContentName(content)
    logger.info("start generate excel file for content :%s" % contshortname)
    excFile = GetScalingResultExcelFile_PerContent(content)
    wb = xlsxwriter.Workbook(excFile)
    shtname = contshortname
    sht = wb.add_worksheet(shtname)

    sht.write(1, 0, 'Content Name')
    sht.write(2, 0, contshortname)
    sht.write(1, 1, 'Scaling Ratio')
    sht.write_column(2, 1, scaleRatios)
    pre_title = ['Width', 'Height', 'DnScaledWidth', 'DnScaledHeight']
    sht.write_row(1, 2, pre_title)
    for dn_algo, up_algo, col in zip(dnScalAlgos, upScalAlgos, ScalQty_WtCols):
        algos = dn_algo + '--' + up_algo
        sht.write(0, col, algos)
        sht.write_row(1, col, QualityList)

    cls, w, h, fr, bitdepth, fmt, totalnum = GetVideoInfo(content, Clips)
    rows = [ScalQty_startRow + i for i in range(len(scaleRatios))]
    continfos = []
    for ratio, row in zip(scaleRatios, rows):
        dw = int(w / ratio)
        dh = int(h / ratio)
        info = [w, h, dw, dh]
        sht.write_row(row, 2, info)
        continfos.append(info)

    charts = []; y_mins = {}; y_maxs = {}
    for qty, x in zip(QualityList, range(len(QualityList))):
        chart_title = 'Scaling Quality - %s' % qty
        xaxis_name = 'scaling ratio'
        chart = CreateChart_Scatter(wb, chart_title, xaxis_name, qty)
        charts.append(chart)
        y_mins[x] = []; y_maxs[x] = []

    for dn_algo, up_algo, col, i in zip(dnScalAlgos, upScalAlgos, ScalQty_WtCols, range(len(dnScalAlgos))):
        qualities = []
        seriname = dn_algo + '--' + up_algo
        for ratio, row, idx in zip(scaleRatios, rows, range(len(scaleRatios))):
            w = continfos[idx][0]
            h = continfos[idx][1]
            dw = continfos[idx][2]
            dh = continfos[idx][3]
            dnScalOut = GetDownScaledOutFile(content, w, h, dw, dh, path_log, dn_algo)  # only care file name
            upScalOut = GetUpScaledOutFile(dnScalOut, dw, dh, w, h, path_log, up_algo)  # only care file name
            qtys = GatherQualityMetrics(upScalOut, path_log)
            sht.write_row(row, col, qtys)
            qualities.append(qtys)
        for x in range(len(QualityList)):
            AddSeriesToChart_Scatter(shtname, rows, col + x, 1, charts[x], seriname, LineColors[i])

        # get min and max of y-axis for a certain dn up scaling algo
        for qty, x in zip(QualityList, range(len(QualityList))):
            qs = [row[x] for row in qualities]
            y_mins[x].append(min(qs))
            y_maxs[x].append(max(qs))

    for qty, x in zip(QualityList, range(len(QualityList))):
        ymin = min(y_mins[x])
        ymax = max(y_maxs[x])
        margin = 0.1  # add 10% on min and max value for y_axis range
        num_precsn = 5 if 'MS-SSIM' in qty else 3
        UpdateChart(charts[x], ymin, ymax, margin, qty, num_precsn)

    startrow = rows[-1] + 2; startcol = 1
    InsertChartsToSheet(sht, startrow, startcol, charts)

    wb.close()
    logger.info("finish export scaling quality results to excel file.")

def GenerateSummarySheet(wb, dnScalAlgos, upScalAlgos, ratio, contentsdict, rows_class, path_log):
    logger.info("start generate summary sheet for ratio %2.2f" % ratio)

    shts = []
    shtname = "Ratio=%1.2f" % ratio
    sht = wb.add_worksheet(shtname)
    shts.append(sht)
    sht.write(1, 0, 'Content Class')
    sht.write(1, 1, 'Content NO.')
    sht.write(1, 2, 'Content Name')
    pre_title = ['Width', 'Height', 'DnScaledWidth', 'DnScaledHeight']
    sht.write_row(1, 3, pre_title)

    for dn_algo, up_algo, col in zip(dnScalAlgos, upScalAlgos, ScalSumQty_WtCols):
        algos = dn_algo + '--' + up_algo
        sht.write(0, col, algos)
        sht.write_row(1, col, QualityList)

    content_infos = {}
    totalnum_contents = 0
    for (clss, contents), row_clss in zip(contentsdict.items(), rows_class):
        sht.write(row_clss, 0, clss)
        totalnum_contents = totalnum_contents + len(contents)
        for content, row_cont in zip(contents, range(len(contents))):
            cl, w, h, fr, bitdepth, fmt, totalnum = GetVideoInfo(content, Clips)
            dw = int(w / ratio)
            dh = int(h / ratio)
            shortcntname = GetShortContentName(content)
            sht.write(row_clss + row_cont, 2, shortcntname)
            infos = [w, h, dw, dh]
            sht.write_row(row_clss + row_cont, 2, infos)
            content_infos[shortcntname] = infos

    charts = []; y_mins = {}; y_maxs = {}
    for qty, x in zip(QualityList, range(len(QualityList))):
        chart_title = '%s of %s' % (qty, shtname)
        xaxis_name = 'contents'
        chart = CreateChart_Scatter(wb, chart_title, xaxis_name, qty)
        charts.append(chart)
        y_mins[x] = []; y_maxs[x] = []

    rows_all = [ScalQty_startRow + i for i in range(totalnum_contents)]
    sht.write_column(ScalQty_startRow, 1, range(totalnum_contents))
    for dn_algo, up_algo, col, i in zip(dnScalAlgos, upScalAlgos, ScalSumQty_WtCols, range(len(dnScalAlgos))):
        qualities = []
        seriname = dn_algo + '--' + up_algo
        for (clss, contents), row_clss in zip(contentsdict.items(), rows_class):
            num_contents = len(contents)
            for content, row_cont in zip(contents, range(num_contents)):
                key = GetShortContentName(content)
                w = content_infos[key][0]
                h = content_infos[key][1]
                dw = content_infos[key][2]
                dh = content_infos[key][3]
                dnScalOut = GetDownScaledOutFile(content, w, h, dw, dh, path_log, dn_algo)  # only care file name
                upScalOut = GetUpScaledOutFile(dnScalOut, dw, dh, w, h, path_log, up_algo)  # only care file name
                qtys = GatherQualityMetrics(upScalOut, path_log)
                sht.write_row(row_clss + row_cont, col, qtys)
                qualities.append(qtys)

        for x in range(len(QualityList)):
            AddSeriesToChart_Scatter(shtname, rows_all, col + x, 1, charts[x], seriname, LineColors[i])

        # get min and max of y-axis for a certain dn up scaling algo
        for qty, x in zip(QualityList, range(len(QualityList))):
            qs = [row[x] for row in qualities]
            y_mins[x].append(min(qs))
            y_maxs[x].append(max(qs))

    for qty, x in zip(QualityList, range(len(QualityList))):
        ymin = min(y_mins[x])
        ymax = max(y_maxs[x])
        margin = 0.1  # add 10% on min and max value for y_axis range
        num_precsn = 5 if 'MS-SSIM' in qty else 3
        UpdateChart(charts[x], ymin, ymax, margin, qty, num_precsn)

    startrow = rows_all[-1] + 2; startcol = 1
    InsertChartsToSheet(sht, startrow, startcol, charts)
    logger.info("finish average sheet for ratio:%2.2f." % ratio)

    return sht

def GenerateAverageSheet(wb, sumsht, dnScalAlgos, upScalAlgos, ratio, contentsdict, rows_class):
    logger.info("start generate average sheet for ratio %2.2f" % ratio)

    rdsht = sumsht
    rdshtname = rdsht.get_name()
    wshtname = "AverageForRatio=%1.2f" % ratio
    wsht = wb.add_worksheet(wshtname)
    wsht.write(1, 0, 'QualityMetrics')
    wsht.write(1, 1, 'Content Class')
    wsht.write(1, 2, 'Content Number')

    interval = 1
    step = len(StatsMetrics) + interval
    startcol = 3
    cols_avg = [startcol + step * i for i in range(len(dnScalAlgos))]
    for col, dn_algo, up_algo in zip(cols_avg, dnScalAlgos, upScalAlgos):
        algos = dn_algo + '--' + up_algo
        wsht.write(0, col, algos)
        wsht.write_row(1, col, StatsMetrics)

    step = len(contentsdict) + 1  # 1 extra row for total of each class
    startrow = 2
    rows_qtymtr = [startrow + step * i for i in range(len(QualityList))]
    for qty, row_qm, y in zip(QualityList, rows_qtymtr, range(len(QualityList))):
        wsht.write(row_qm, 0, qty)

        #charts = []
        #titlename = 'Quality Statistics %s' % qty
        #chart = CreateChart_Line(wb, titlename, qty)
        #charts.append(chart)

        totalnum_contents = 0
        for (cls, contents), idx, rdrow_cls in zip(contentsdict.items(), range(len(contentsdict)), rows_class):
            wsht.write(row_qm + idx, 1, cls)
            num_content = len(contents)
            totalnum_contents = totalnum_contents + num_content
            wsht.write(row_qm + idx, 2, num_content)
            for rdcol, wtcol in zip(ScalSumQty_WtCols, cols_avg):
                startcell = xlrd.cellname(rdrow_cls, rdcol + y)
                endcell = xlrd.cellname(rdrow_cls + num_content - 1, rdcol + y)
                formula = '=MAX(\'%s\'!%s:%s)' % (rdshtname, startcell, endcell)
                wsht.write(row_qm + idx, wtcol, formula)

                formula = '=SUM(\'%s\'!%s:%s)/%d' % (rdshtname, startcell, endcell, num_content)
                wsht.write(row_qm + idx, wtcol + 2, formula)

                formula = '=MIN(\'%s\'!%s:%s)' % (rdshtname, startcell, endcell)
                wsht.write(row_qm + idx, wtcol + 4, formula)

                avgcell = xlrd.cellname(row_qm + idx, wtcol + 2)
                formula = '= %s + _xlfn.STDEV.P(\'%s\'!%s:%s)' % (avgcell, rdshtname, startcell, endcell)
                wsht.write(row_qm + idx, wtcol + 1, formula)

                formula = '= %s - _xlfn.STDEV.P(\'%s\'!%s:%s)' % (avgcell, rdshtname, startcell, endcell)
                wsht.write(row_qm + idx, wtcol + 3, formula)

        #write total contents statistics
        wtrow = row_qm + len(contentsdict)
        wsht.write(wtrow, 1, 'Total')
        wsht.write(wtrow, 2, totalnum_contents)
        for rdcol, wtcol in zip(ScalSumQty_WtCols, cols_avg):
            startcell = xlrd.cellname(ScalQty_startRow, rdcol + y)
            endcell = xlrd.cellname(ScalQty_startRow + totalnum_contents - 1, rdcol + y)
            formula = '=MAX(\'%s\'!%s:%s)' % (rdshtname, startcell, endcell)
            wsht.write(wtrow, wtcol, formula)

            formula = '=SUM(\'%s\'!%s:%s)/%d' % (rdshtname, startcell, endcell, totalnum_contents)
            wsht.write(wtrow, wtcol + 2, formula)

            formula = '=MIN(\'%s\'!%s:%s)' % (rdshtname, startcell, endcell)
            wsht.write(wtrow, wtcol + 4, formula)

            avgcell = xlrd.cellname(wtrow, wtcol + 2)
            formula = '= %s + _xlfn.STDEV.P(\'%s\'!%s:%s)' % (avgcell, rdshtname, startcell, endcell)
            wsht.write(wtrow, wtcol + 1, formula)

            formula = '= %s - _xlfn.STDEV.P(\'%s\'!%s:%s)' % (avgcell, rdshtname, startcell, endcell)
            wsht.write(wtrow, wtcol + 3, formula)

    logger.info("finish average sheet for ratio:%2.2f." % ratio)


def SaveScalingResultsToExcel(dnScalAlgos, upScalAlgos, path_log):
    logger.info("start saving scaling quality results to excel files.......")
    if not os.path.exists(Path_ScalingResults):
        os.makedirs(Path_ScalingResults)

    contents = GetContents(ContentPath, Clips)
    for content in contents:
        GeneratePerContentExcelFile(dnScalAlgos, upScalAlgos, content, DnScaleRatio, path_log)

    scaleRatios = DnScaleRatio
    scaleRatios.remove(1.0)
    logger.info("start generating scaling quality summary excel file.......")
    sumexcFile = GetScalingResultExcelFile(len(scaleRatios), len(DnScaleRatio))
    wb = xlsxwriter.Workbook(sumexcFile)
    # to generate rows number of starting of each class: rows_class
    contentsdict, rows_class = CalcRowsClassAndContentDict(ScalQty_startRow, ContentPath, Clips)

    for ratio in scaleRatios:
        sht = GenerateSummarySheet(wb, dnScalAlgos, upScalAlgos, ratio, contentsdict, rows_class, path_log)
        GenerateAverageSheet(wb, sht, dnScalAlgos, upScalAlgos, ratio, contentsdict, rows_class)

    wb.close()
    logger.info("finish export saving scaling quality results to excel files.......")


