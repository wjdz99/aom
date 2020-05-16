This testing framework is initial developed by Intel for convex hull study. It is contributed to Alliance of
Open Media for convex hull study and future AV2 Common Test Condition (CTC) testing.

For questions and technical support, please contact ryan.lei@intel.com or maggie.sun@intel.com

This readme file provides a simple introduction of the framework and steps to use it for different studies.

Prerequisites
1. Inside ./bin folder, make sure you have excutables for all external tools, such as vmafossexec, ffmpeg,
   HDRMetrics, encoder, decoder, etc. If you need to run on Linux, you need to get the linux executables
   for those tools.
2. VMAF requires a model file to run. currently, vmaf_v0.6.1 model files (vmaf_v0.6.1.pkl.model, vmaf_v0.6.1.pkl)
   are used and stored under the ./bin folder.
3. HDRMetrics requires a template config file (HDRMetricsYUV420_template.cfg). It is also located under the
   ./bin folder. Test framework will generate individual config file based on this template.
4. VBA macro binary file (vbaProject-AV2.bin) for bdrate calculation is also stored under ./bin folder, which
   is needed for BDRATE calculation inside the excel file.
5. Test framework is implemented in python 3. It requires few extra python packages, such as xlrd, xlsxwriter,
   argparse, etc.

Things you need to update inside ./src/Config.py to configure the test.
1. Update the Clips table following the existing example to provide the list of video sequences that you want
   to test. The first value is the yuv file name without .yuv suffix. All test video sequences should be located
   in the folder specified by ContentPath. Test framework will extract the string before the first "_" as the
   short name to identify the test sequence and use that to generate file name of other intermediate files, so
   please make sure the short name of the test sequence is unique. Currently only 8 bit yuv sequences in yuv420p
   format are supported.

2. Update FrameNum to specify number of frames you want to process.

3. Update DnScaleRatio list to provide the list of downscaling ratios.

4. Update DnScalingAlgos and UpScalingAlgos list to specify downscaling/upscaling filter types that you want to
   test. The name of the filter types are that supported by ffmpeg. The size of these two lists have to be the same.
   Filter types for downscaling and upscaling can be different.

5. Update QPs list to specify the list of QPs you want to test for different codecs. QPs have to be in the valid QP
   range for different coding standards. For example, [0, 51] for hevc and [0, 63] for AV1

6. Update QualityEvalMethod to specify tools you want to use to calculate quality metrics. Currently, only VMAF and
   HDRTools are supported. You also need to update QualityList accordingly to specify what quality metrics you want
   to calculate. Test framework will filter out other metrics and only report the metrics that are specified.

7. Update LogCmdOnly flag. When it is set to True, the test framework will only capture all process command sequences
   into a log file (under ./test/testLogs folder with time stamp) without actually running it on your system. This
   feature is to support the use case in which actual processing tasks need to distributed on to a server cluster and
   executed in parallel.

8. Update SMOKE_TEST flag, It is used for sanity check purpose, in which only few frames are processed to reduce the
   test time.

Command lines for running the tests:
below is the full command line options in help message:
  -h, --help                   show this help message and exit
  -f, --function               function to run: clean, scaling, sumscaling, encode, convexhull, summary
  -k, --KeepUpscaleOutput      [0|1] in function clean, if keep upscaled yuv files
  -s, --SaveMemory             [0|1] save memory mode will delete most files in intermediate steps and keeps only
                               necessary ones for RD calculation
  -l, --LoggingLevel           logging level: 0:No Logging, 1: Critical, 2: Error, 3:Warning, 4: Info, 5: Debug
  -c, --CodecName              CodecName: av1 or hevc
  -m, --EncodeMethod           EncodeMethod: ffmpeg, aom, svt
  -p, --EncodePreset           EncodePreset: medium or slow for ffmpeg, 0,1,2... for aom and svt

Sample command for typical operations:
1.  python ConvexHullTest.py -f clean
    This command will clean up all intermediate files under ./test folder

2.  python ConvexHullTest.py -f scaling
    This command will run the standalone downscaling and upscaling tests. Downscaled YUV files are stored under
    ./test/downscaledYUVs folder. Upscaled YUV files are stored under ./test/upscaledYUVs folder. Quality metrics
    log files are stored under ./test/qualityLogs folder. Other processing logs and command logs are stored under
    ./test/testLogs folder. In case HDRMetric is used, individual config files are generated and stored under
    ./test/configFiles folder. All intermediate file names indicate the input, output resolution and filter types
    that are used.

3.  python ConvexHullTest.py -f sumscaling
    This command will summarize the quality metrics for the scaling test into excel files under ./analysis/scalingresult
    folder. There are excel files for each individual test sequence and also excel file that summarizes quality result
    for all test sequences based on classes.

4.  python ConvexHullTest.py -f encode -c hevc -m ffmpeg -p medium
    This command will run the encoding test. It actually contains encode, decode, upscale, quality metrics steps. Different
    encoding method and codec name can be specified. For a particular encoder, encoding preset should also be specified.
    Encoded bitstreams are stored under ./test/bitstreams folder. Decoded YUV files are stored under ./test/decodedYUVs
    folder. Decoded and then upscaled YUV files are stored under ./test/decUpscaledYUVs folder. Quality logs are stored
    under ./test/qualityLogs folder.

5.  python ConvexHullTest.py -f convexhull -c hevc -m ffmpeg -p medium
    This command will summarize the per sequence quality result based on different scaling ratios and scaling filter types.
    Please make sure the same encoding method/codec name/preset are used as the previous steps. Output excel files are
    stored under ./analysis/rdresult folder. For each sequence and downscaling/upscaling filter type, an excel sheet is
    generated that contains the bitrate and quality metric for different downscaled resolutions. Rate distortion curve
    for all quality metrics will be drawn in a scatter plot. Convex hull will be calculated and draw on top of the rate
    distortion curve. BDRATE between different resolutions will be calculated.

6.  python ConvexHullTest.py -f summary -c hevc -m ffmpeg -p medium
    This command will summarize the quality metrics across all test sequences into an excel file stored under the
    ./analysis/summary folder. Average result based on content classes will be calculated.

