import numpy as np

x = np.genfromtxt("data/lowres_64f_target400_data.txt", delimiter="\t")

prev_update = 0
split_list_indices = list()
for i, val in enumerate(x):
    if prev_update != val[2]:
        split_list_indices.append(i)
        prev_update = val[2]

split = np.split(x, split_list_indices)

for array in split:
    A, B, update = np.hsplit(array, 3)
    x, residuals, rank, s = np.linalg.lstsq(A, B, rcond=None)
    print(update[0], x)
