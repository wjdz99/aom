import numpy as np

data = np.genfromtxt("data/lowres_64f_target150_data.txt", delimiter="\t")

prev_update = 0
split_list_indices = list()
for i, val in enumerate(data):
    if prev_update != val[3]:
        split_list_indices.append(i)
        prev_update = val[3]

split = np.split(data, split_list_indices)

for array in split:
    A, mv, B, update = np.hsplit(array, 4)
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    new_A = np.concatenate((A, mv), axis=1)
    new_A_inv = np.linalg.pinv(new_A)
    new_x = np.matmul(new_A_inv, B)
    print("update type:", update[0][0])
    print("least squares solution:", x[0][0])
    print("pinv solution:", new_x[0][0], new_x[1][0])
    print()
