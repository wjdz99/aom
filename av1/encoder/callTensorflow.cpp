﻿#include "/usr/include/python3.6m/Python.h"
#include "stdio.h"
#include <sstream>

#include <limits.h>
#include <math.h>

#include "config/aom_config.h"
#include "config/aom_dsp_rtcd.h"
#include "config/aom_scale_rtcd.h"
#include "config/av1_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/aom_filter.h"
#if CONFIG_DENOISE
#include "aom_dsp/grain_table.h"
#include "aom_dsp/noise_util.h"
#include "aom_dsp/noise_model.h"
#endif
#include "aom_dsp/psnr.h"
#if CONFIG_INTERNAL_STATS
#include "aom_dsp/ssim.h"
#endif
#include "aom_ports/aom_timer.h"
#include "aom_ports/mem.h"
#include "aom_ports/system_state.h"
#include "aom_scale/aom_scale.h"
#if CONFIG_BITSTREAM_DEBUG || CONFIG_MISMATCH_DEBUG
#include "aom_util/debug_util.h"
#endif  // CONFIG_BITSTREAM_DEBUG || CONFIG_MISMATCH_DEBUG

#include "av1/common/alloccommon.h"
#include "av1/common/cdef.h"
#include "av1/common/filter.h"
#include "av1/common/idct.h"
#include "av1/common/reconinter.h"
#include "av1/common/reconintra.h"
#include "av1/common/resize.h"

#include "av1/encoder/aq_complexity.h"
#include "av1/encoder/aq_cyclicrefresh.h"
#include "av1/encoder/aq_variance.h"
#include "av1/encoder/bitstream.h"
#include "av1/encoder/context_tree.h"
#include "av1/encoder/encodeframe.h"
#include "av1/encoder/encodemv.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/encodetxb.h"
#include "av1/encoder/ethread.h"
#include "av1/encoder/firstpass.h"
#include "av1/encoder/grain_test_vectors.h"
#include "av1/encoder/hash_motion.h"
#include "av1/encoder/mbgraph.h"
#include "av1/encoder/picklpf.h"
#include "av1/encoder/pickrst.h"
#include "av1/encoder/random.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/speed_features.h"
#include "av1/encoder/temporal_filter.h"


using namespace std;

/*
  author：cgy
  date：2018/9/22
  param：
  ppp： 指向 encoder.c文件中cm->frame_to_show->y_buffer的地址，即重建图像的首地址；
  height：序列的高；
  width：序列的宽；

*/

uint8_t** callTensorflow(uint8_t* ppp, int height, int width, int stride, FRAME_TYPE frame_type) {

  Py_SetPath(L"/usr/local/google/home/logangw/aom/av1/encoder:"
             "/usr/lib:"
             "/usr/lib/python3.6:"
             "/usr/lib/python3.6/site-packages:"
             "/usr/lib/python3.6/lib-dynload");

  PyObject* pModule = NULL;
  PyObject* pFuncI = NULL;
  PyObject* pFuncB = NULL;
  PyObject* pArgs = NULL;

  // 初始化python环境
  Py_Initialize();

  if (!Py_IsInitialized()) {
    printf("Python init failed!\n");
    return NULL;
  }

  // import python
  pModule = PyImport_ImportModule("TEST");

  if (!pModule) {
    printf("don't load Pmodule\n");
    Py_Finalize();
    return NULL;
  }
  // 获得TensorFlow函数指针
  pFuncI = PyObject_GetAttrString(pModule, "entranceI");
  if (!pFuncI) {
    printf("don't get I function!");
    Py_Finalize();
    return NULL;
  }

  pFuncB = PyObject_GetAttrString(pModule, "entranceB");
  if (!pFuncB) {
    printf("don't get B function!");
    Py_Finalize();
    return NULL;
  }

  PyObject* list = PyList_New(height);
  pArgs = PyTuple_New(1);                 //以元组方式传参
  PyObject** lists = new PyObject*[height];

  //将图像缓冲区的数据读到列表中
  for (int i = 0; i < height; i++)
  {
    lists[i] = PyList_New(0);
    for (int j = 0; j < width; j++)
    {
      PyList_Append(lists[i], Py_BuildValue("i", *(ppp + j)));//转化为python对象
    }
    PyList_SetItem(list, i, lists[i]);
    ppp += stride;
  }
  PyTuple_SetItem(pArgs, 0, list);    //"list" is the input image

  PyObject *presult = NULL;

  if (frame_type == KEY_FRAME){
    presult = PyEval_CallObject(pFuncI, pArgs);//将pArgs参数传递到Python中
  }
  else{
    presult = PyEval_CallObject(pFuncB, pArgs);//将pArgs参数传递到Python中
  }

  //需要定义一个二维数组；
  uint8_t **rePic = new uint8_t*[height];
  for (int i = 0; i < height; i++){
    rePic[i] = new uint8_t[width];
  }
  uint8_t s;

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++)  {
      PyArg_Parse(PyList_GetItem(PyList_GetItem(presult, i), j), "B", &s);
      rePic[i][j] = s;
    }
  }

  return rePic;
}


uint16_t** callTensorflow_hbd(uint16_t* ppp, int height, int width, int stride, FRAME_TYPE frame_type){

  Py_SetPath(L"/usr/local/google/home/logangw/aom/av1/encoder:"
             "/usr/lib:"
             "/usr/lib/python3.6:"
             "/usr/lib/python3.6/site-packages:"
             "/usr/lib/python3.6/lib-dynload");

  PyObject* pModule = NULL;
  PyObject* pFuncI = NULL;
  PyObject* pFuncB = NULL;
  PyObject* pArgs = NULL;

  // 初始化python环境
  Py_Initialize();

  if (!Py_IsInitialized()) {
    printf("Python init failed!\n");
    return NULL;
  }

  pModule = PyImport_ImportModule("TEST");

  if (!pModule) {
    printf("don't load Pmodule\n");
    Py_Finalize();
    return NULL;
  }

  // 获得TensorFlow函数指针
  pFuncI = PyObject_GetAttrString(pModule, "entranceI");
  if (!pFuncI) {
    printf("don't get I function!");
    Py_Finalize();
    return NULL;
  }

  pFuncB = PyObject_GetAttrString(pModule, "entranceB");
  if (!pFuncB) {
    printf("don't get B function!");
    Py_Finalize();
    return NULL;
  }

  PyObject* list = PyList_New(height);
  pArgs = PyTuple_New(1);                 //以元组方式传参
  PyObject** lists = new PyObject*[height];

  //Read the data from y buffer into the list
  for (int i = 0; i < height; i++)
  {
    lists[i] = PyList_New(0);
    for (int j = 0; j < width; j++)
    {
      PyList_Append(lists[i], Py_BuildValue("i", *(ppp + j)));//Convert to Python objects
    }
    PyList_SetItem(list, i, lists[i]);
    ppp += stride;
  }
  PyTuple_SetItem(pArgs, 0, list);    //将列表赋给参数

  PyObject *presult = NULL;

  if (frame_type == KEY_FRAME){
    presult = PyEval_CallObject(pFuncI, pArgs);//将pArgs参数传递到Python中
  }
  else{
    presult = PyEval_CallObject(pFuncB, pArgs);//将pArgs参数传递到Python中
  }

  //需要定义一个二维数组；
  uint16_t **rePic = new uint16_t*[height];
  for (int i = 0; i < height; i++){
    rePic[i] = new uint16_t[width];
  }
  uint16_t s;

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      PyArg_Parse(PyList_GetItem(PyList_GetItem(presult, i), j), "H", &s);
      rePic[i][j] = s;
    }
  }

  //关闭python解释器
  return rePic;
}


uint8_t** blockCallTensorflow(uint8_t* ppp, int cur_buf_height, int cur_buf_width, int stride, FRAME_TYPE frame_type) {

  Py_SetPath(L"/usr/local/google/home/logangw/aom/av1/encoder:"
             "/usr/lib:"
             "/usr/lib/python3.6:"
             "/usr/lib/python3.6/site-packages:"
             "/usr/lib/python3.6/lib-dynload");

  PyObject * pModule = NULL;
  PyObject * pFuncI = NULL;
  PyObject * pFuncB = NULL;
  PyObject * pArgs = NULL;

  Py_Initialize();

  if (!Py_IsInitialized()) {
    printf("Python init failed!\n");
    return NULL;
  }

  pModule = PyImport_ImportModule("TEST");

  if (!pModule) {
    printf("don't load Pmodule\n");
    Py_Finalize();
    return NULL;
  }

  pFuncI = PyObject_GetAttrString(pModule, "entranceI");
  if (!pFuncI) {
    printf("don't get I function!");
    Py_Finalize();
    return NULL;
  }

  pFuncB = PyObject_GetAttrString(pModule, "entranceB");
  if (!pFuncB) {
    printf("don't get B function!");
    Py_Finalize();
    return NULL;
  }

  PyObject* list = PyList_New(cur_buf_height);
  pArgs = PyTuple_New(1);
  PyObject** lists = new PyObject*[cur_buf_height];

  for (int i = 0; i < cur_buf_height; i++)
  {
    lists[i] = PyList_New(0);
    for (int j = 0; j < cur_buf_width; j++)
    {
      PyList_Append(lists[i], Py_BuildValue("i", *(ppp + j)));
    }
    PyList_SetItem(list, i, lists[i]);
    ppp += stride;
  }

  PyTuple_SetItem(pArgs, 0, list);
  PyObject *presult = NULL;
  if (frame_type == KEY_FRAME){
    presult = PyEval_CallObject(pFuncI, pArgs);
  } else {
    presult = PyEval_CallObject(pFuncB, pArgs);
  }

  uint8_t **rePic = new uint8_t*[cur_buf_height];
  for (int i = 0; i < cur_buf_height; i++) {
    rePic[i] = new uint8_t[cur_buf_width];
  }
  uint8_t s;

  for (int i = 0; i < cur_buf_height; i++)
  {
    for (int j = 0; j < cur_buf_width; j++)
    {
      PyArg_Parse(PyList_GetItem(PyList_GetItem(presult, i), j), "i", &s);
      rePic[i][j] = s;
    }
  }
  return rePic;
}


uint16_t** blockCallTensorflow_hbd(uint16_t* ppp, int cur_buf_height, int cur_buf_width, int stride, FRAME_TYPE frame_type) {

  Py_SetPath(L"/usr/local/google/home/logangw/aom/av1/encoder:"
             "/usr/lib:"
             "/usr/lib/python3.6:"
             "/usr/lib/python3.6/site-packages:"
             "/usr/lib/python3.6/lib-dynload");

  PyObject* pModule = NULL;
  PyObject* pFuncI = NULL;
  PyObject* pFuncB = NULL;
  PyObject* pArgs = NULL;

  Py_Initialize();

  if (!Py_IsInitialized()) {
    printf("Python init failed!\n");
    return NULL;
  }

  pModule = PyImport_ImportModule("TEST");

  if (!pModule) {
    printf("don't load Pmodule\n");
    Py_Finalize();
    return NULL;
  }

  pFuncI = PyObject_GetAttrString(pModule, "entranceI");
  if (!pFuncI) {
    printf("don't get I function!");
    Py_Finalize();
    return NULL;
  }
  pFuncB = PyObject_GetAttrString(pModule, "entranceB");
  if (!pFuncB) {
    printf("don't get B function!");
    Py_Finalize();
    return NULL;
  }
  PyObject* list = PyList_New(cur_buf_height);
  pArgs = PyTuple_New(1);
  PyObject** lists = new PyObject*[cur_buf_height];

  for (int i = 0; i < cur_buf_height; i++)
  {
    lists[i] = PyList_New(0);
    for (int j = 0; j < cur_buf_width; j++)
    {
      PyList_Append(lists[i], Py_BuildValue("i", *(ppp + j)));
    }
    PyList_SetItem(list, i, lists[i]);
    ppp += stride;
  }
  PyTuple_SetItem(pArgs, 0, list);
  PyObject *presult = NULL;
  if (frame_type == KEY_FRAME){
    presult = PyEval_CallObject(pFuncI, pArgs);
  }
  else{
    presult = PyEval_CallObject(pFuncB, pArgs);
  }

  uint16_t **rePic = new uint16_t*[cur_buf_height];
  for (int i = 0; i < cur_buf_height; i++) {
    rePic[i] = new uint16_t[cur_buf_width];
  }
  uint16_t s;

  for (int i = 0; i < cur_buf_height; i++) {
    for (int j = 0; j < cur_buf_width; j++) {
      PyArg_Parse(PyList_GetItem(PyList_GetItem(presult, i), j), "H", &s);
      rePic[i][j] = s;
    }
  }
  return rePic;
}
