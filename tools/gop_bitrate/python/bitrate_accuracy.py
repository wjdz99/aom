import numpy as np

# Uses least squares regression to find the solution
# when there is one unknown variable.
def lstsq_solution(A, B):
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    return x[0][0]

# Uses the pseudoinverse matrix to find the solution
# when there are two unknown variables.
def pinv_solution(A, mv, B):
    new_A = np.concatenate((A, mv), axis=1)
    new_A_inv = np.linalg.pinv(new_A)
    new_x = np.matmul(new_A_inv, B)
    print("pinv solution:", new_x[0][0], new_x[1][0])
    return (new_x[0][0], new_x[1][0])

# Calculates the least squares error between A and B
# using coefficients in X.
def lstsq_error(A, B, x):
    error = 0
    n = 0
    for i, a in enumerate(A):
        a = a[0]
        b = B[i][0]
        if b == 0:
            continue
        n += 1
        error += (b - x*a)**2
    if n == 0:
        return None
    error /= n
    return error

# Calculates the percentage error between A and B
# using coefficients in X.
def percent_error(A, B, x):
    error = 0
    n = 0
    for i, a in enumerate(A):
        a = a[0]
        b = B[i][0]
        if b == 0:
            continue
        n += 1
        error += (abs(x*a-b)/b)*100
    if n == 0:
        print("hi")
        return None
    error /= n
    return error

# Traverses the data and prints out one value for
# each update type.
def print_solutions(file_path):
    data = np.genfromtxt(file_path, delimiter="\t")
    prev_update = 0
    split_list_indices = list()
    for i, val in enumerate(data):
        if prev_update != val[3]:
            split_list_indices.append(i)
            prev_update = val[3]
    split = np.split(data, split_list_indices)
    for array in split:
        A, mv, B, update = np.hsplit(array, 4)
        z = np.where(B == 0)[0]
        r_e = np.delete(A, z, axis=0)
        r_m = np.delete(mv, z, axis=0)
        r_f = np.delete(B, z, axis=0)
        A = r_e
        mv = r_m
        B = r_f
        all_zeros = not A.any()
        if all_zeros:
            continue
        print("update type:", update[0][0])
        x = lstsq_solution(A, B)
        print("x:", x)
        training_error = lstsq_error(A, B, x)
        baseline_error = lstsq_error(A, B, 1)
        percent_training_error = percent_error(A, B, x)
        percent_baseline_error = percent_error(A, B, 1)
        print("lstsq training error:", training_error, "lstsq baseline error:", baseline_error)
        percent_lstsq_error = (abs(training_error - baseline_error)/baseline_error)*100
        print("percent training error:", percent_training_error, "percent baseline error:", percent_baseline_error)
        print()

if __name__ == "__main__":
    print_solutions("data2/lowres_17f_target_gt400_data.txt")
