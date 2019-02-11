import sys
import numpy as np
import tensorflow as tf
import os, time
from VDSR15 import model15
from VDSR20 import model20
from VDSR25 import model25
from VDSR30 import model30
from UTILS import *

if not hasattr(sys, 'argv'):
  sys.argv = ['']

I_MODEL_PATH = r"/usr/local/google/home/logangw/aom/models/I_model"  # epoch 364
B_MODEL_PATH = r"/usr/local/google/home/logangw/aom/models/B_model"  # epoch 252


def prepare_test_data(fileOrDir):
  original_ycbcr = []
  gt_y = []
  fileName_list = []
  imgCbCr = 0
  fileName_list.append(fileOrDir)
  imgY = np.reshape(fileOrDir,(1, len(fileOrDir), len(fileOrDir[0]), 1))
  imgY = normalize(imgY)

  original_ycbcr.append([imgY, imgCbCr])
  return original_ycbcr, gt_y, fileName_list


def test_all_ckpt(modelPath, fileOrDir, is_key_frame_model):
  tf.reset_default_graph()
  tf.logging.warning(modelPath)

  tem = [f for f in os.listdir(modelPath) if 'data' in f]
  ckptFiles = sorted([r.split('.data')[0] for r in tem])

  config = tf.ConfigProto()
  config.gpu_options.allow_growth = True
  with tf.Session(config=config) as sess:
    input_tensor = tf.placeholder(tf.float32, shape=(1, None, None, 1))
    if is_key_frame_model:
      shared_model = tf.make_template('shared_model', model25)
    else:
      shared_model = tf.make_template('shared_model', model25)

    output_tensor, weights = shared_model(input_tensor)
    output_tensor = tf.clip_by_value(output_tensor, 0., 1.)
    output_tensor = output_tensor * 255

    sess.run(tf.global_variables_initializer())

    original_ycbcr, gt_y, fileName_list = prepare_test_data(fileOrDir)

    for ckpt in ckptFiles:
      tf.logging.warning(ckpt)
      epoch = int(ckpt.split('_')[-1].split('.')[0])

      if is_key_frame_model:
        if epoch != 364:
          continue
      else:
        if epoch != 252:
          continue

      saver = tf.train.Saver(tf.global_variables())
      saver.restore(sess, os.path.join(modelPath, ckpt))
      total_imgs = len(fileName_list)
      for i in range(total_imgs):
        imgY = original_ycbcr[i][0]
        out = sess.run(output_tensor, feed_dict={input_tensor: imgY})
        out = np.reshape(out, (out.shape[1], out.shape[2]))
        out = np.around(out)
        out = out.astype('int')
        out = out.tolist()

        return out


def entranceI(inp):
  tf.logging.warning("python, in I")
  i = test_all_ckpt(I_MODEL_PATH, inp, True)
  return i


def entranceB(inp):
  tf.logging.warning("python, in B")
  b = test_all_ckpt(B_MODEL_PATH, inp, False)
  return b
