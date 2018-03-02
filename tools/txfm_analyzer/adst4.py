from math import *
import numpy as np

def gen_adst4():
  A = np.zeros((4, 4))
  for j in range(1, 5):
    for i in range(1, 5):
      A[j-1, i-1] = sin((2*j - 1) * i * pi/9) * 2 / 3 * (2**0.5)
  return A

def gen_sinpi():
  for bit in range(10, 17):
    for k in range(5):
      sinpi = (2**bit) * (2**0.5) * sin(k * pi / 9) * 2 / 3
      sinpi = int(round(sinpi))
      print sinpi, ",",
    print ""

def sin_sc(t):
  sc = 2. / 3 * (2**0.5)
  return sc*sin(t)

def gen_fadst4_matrix():
  M_ls = []

  # stage 0
  X = np.identity(4)
  M_ls.append(X)

  # stage 1
  A = [[sin_sc(pi/9)  , 0          , 0          , 0],
       [sin_sc(4*pi/9), 0          , 0          , 0],
       [0          , sin_sc(2*pi/9), 0          , 0],
       [0          , sin_sc(pi/9)  , 0          , 0],
       [0          , 0          , sin_sc(3*pi/9), 0],
       [0          , 0          , 0          , sin_sc(4*pi/9)],
       [0          , 0          , 0          , sin_sc(2*pi/9)],
       [1          , 1          , 0          , 0],
       [0          , 0          , 0          , -1]];
  A = np.array(A)
  M_ls.append(A)

  # stage 2
  B = np.zeros((8, 9))
  for i in range(8):
    B[i, i] = 1;
  B[7, 8] = 1
  M_ls.append(B)

  # stage 3
  C = np.zeros((6, 8))

  C[0, 5] = 1

  C[1, 0] = 1
  C[1, 2] = 1

  C[2, 7] = sin_sc(3*pi/9)

  C[3, 6] = 1

  C[4, 1] = 1
  C[4, 3] = -1

  C[5, 4] = 1

  M_ls.append(C)

  # stage 4
  D = np.zeros((4, 6))
  D[0, 0] = 1
  D[0, 1] = 1
  D[1, 2] = 1
  D[2, 3] = 1
  D[2, 4] = 1
  D[3, 5] = 1

  M_ls.append(D)

  # stage 5
  E = np.zeros((5, 4))
  E[0, 0] = 1
  E[0, 3] = 1
  E[1, 1] = 1
  E[2, 2] = 1
  E[2, 3] = -1
  E[3, 0] = -1
  E[3, 2] = 1
  E[4, 3] = 1

  M_ls.append(E)

  # stage 6
  F = np.zeros((4, 5))
  F[0, 0] = 1
  F[1, 1] = 1
  F[2, 2] = 1
  F[3, 3] = 1
  F[3, 4] = 1

  M_ls.append(F)
  X = np.identity(4)
  for M in M_ls:
    X = np.dot(M, X)
  return X, M_ls

def gen_iadst4_matrix():
  M_ls = []

  # stage 0
  X = np.identity(4)
  M_ls.append(X)

  # stage 1
  A = [[sin_sc(pi/9)  , 0          , 0          , 0], #sinpi[1] * x0
       [sin_sc(2*pi/9), 0          , 0          , 0], #sinpi[2] * x0
       [0          , sin_sc(3*pi/9), 0          , 0], #sinpi[3] * x1
       [0          , 0          , sin_sc(4*pi/9), 0], #sinpi[4] * x2
       [0          , 0          , sin_sc(pi/9)  , 0], #sinpi[1] * x2
       [0          , 0          , 0          , sin_sc(2*pi/9)], #sinpi[2] * x3
       [0          , 0          , 0          , sin_sc(4*pi/9)], #sinpi[3] * x3
       [1          , 0          , -1         , 0], #x0 - x2
       [0          , 0          , 0          , 1]]; #x3
  A = np.array(A)
  M_ls.append(A)

  # stage 2
  B = np.zeros((8, 9))
  for i in range(8):
    B[i, i] = 1;
  B[7, 8] = 1 # B[7, 7] = 1 and B[7, 8] = 1 -> (x0 - x2) + x3
  M_ls.append(B)

  # stage 3
  C = np.zeros((6, 8))

  C[0, 5] = 1 # s5

  C[1, 0] = 1 # s0 = s0 + s3
  C[1, 3] = 1 # s0 = s0 + s3

  C[2, 6] = -1 # -s6

  C[3, 1] = 1 # s1 = s1 - s4
  C[3, 4] = -1 # s1 = s1 - s4

  C[4, 7] = sin_sc(3*pi/9) # s2 = sinpi[3] * s7

  C[5, 2] = 1 # s3 = s2


  M_ls.append(C)

  # stage 4
  D = np.zeros((4, 6))
  D[0, 0] = 1 # s0 = s0 + s5
  D[0, 1] = 1 # s0 = s0 + s5
  D[1, 2] = 1 # s1 = s1 - s6
  D[1, 3] = 1 # s1 = s1 - s6
  D[2, 4] = 1 # s2 = -s6
  D[3, 5] = 1 # s3

  M_ls.append(D)

  # stage 5
  E = np.zeros((5, 4))
  E[0, 0] = 1 #x0 = s0 + s3
  E[0, 3] = 1 #x0 = s0 + s3
  E[1, 1] = 1 #x1 = s1 + s3
  E[1, 3] = 1 #x1 = s1 + s3
  E[2, 2] = 1 #x2 = s2
  E[3, 0] = 1 #x3 = s0 + s1
  E[3, 1] = 1 #x3 = s0 + s1
  E[4, 3] = -1 # -s3

  M_ls.append(E)

  # stage 6
  F = np.zeros((4, 5))
  F[0, 0] = 1 # x0
  F[1, 1] = 1 # x1
  F[2, 2] = 1 # x2
  F[3, 3] = 1 # x3 - s3
  F[3, 4] = 1 # x3 - s3

  M_ls.append(F)
  X = np.identity(4)
  for M in M_ls:
    X = np.dot(M, X)
  return X, M_ls

if __name__ == "__main__":
  #A = gen_adst4()
  #print A
  X_fadst, M_fadst = gen_fadst4_matrix()
  X_iadst, M_iadst = gen_iadst4_matrix()
  M_ls = M_fadst + M_iadst

  X = np.identity(4)
  print "fadst stage range"
  for i, M in enumerate(M_ls[0:7]):
    X = np.dot(M, X)
    m = max(np.sum(abs(X), axis = 1))
    print "stage", i, "bit range", log(m, 2)

  print "iadst stage range"
  for i, M in enumerate(M_ls[7:14]):
    X = np.dot(M, X)
    m = max(np.sum(abs(X), axis = 1))
    print "stage", i, "bit range", log(m, 2)
