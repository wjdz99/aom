import numpy as np

# Removes indices from both arrays where zeros occur.
def remove_zero_entries_sync(A, B):
    z = np.where(B == 0)[0]
    A = np.delete(A, z, axis=0)
    B = np.delete(B, z, axis=0)
    z = np.where(A == 0)[0]
    A = np.delete(A, z, axis=0)
    B = np.delete(B, z, axis=0)
    return (A, B)

# Finds the coefficient to multiply A by to minimize
# the percentage error between A and B.
def minimize_percentage_error_model_a(A, B):
    z = np.where(B == 0)[0]
    A = np.delete(A, z, axis=0)
    B = np.delete(B, z, axis=0)
    z = np.where(A == 0)[0]
    A = np.delete(A, z, axis=0)
    B = np.delete(B, z, axis=0)

    R = np.divide(A, B)
    num = 0
    den = 0
    for r_i in R:
        num += r_i
        den += r_i**2
    if den == 0:
        x = 0
    else:
        x = (num / den)[0]
    return x

def minimize_percentage_error_model_b(r_e, r_m, r_f):
    z = np.where(r_f == 0)[0]
    r_e = np.delete(r_e, z, axis=0)
    r_m = np.delete(r_m, z, axis=0)
    r_f = np.delete(r_f, z, axis=0)

    r_ef = np.divide(r_e, r_f)
    r_mf = np.divide(r_m, r_f)
    sum_ef = np.sum(r_ef)
    sum_ef_sq = np.sum(np.square(r_ef))
    sum_mf = np.sum(r_mf)
    sum_mf_sq = np.sum(np.square(r_mf))
    sum_ef_mf = np.sum(np.multiply(r_ef, r_mf))

    # Sometimes all the estimated MVs are 0.
    # In this case we can revert to model A.
    if sum_mf_sq == 0 or sum_ef_sq == 0:
        return (minimize_percentage_error_model_a(r_e, r_f), 1)

    # Set up and solve the matrix equation
    A = np.array([[1, (sum_ef_mf / sum_ef_sq)],[(sum_ef_mf / sum_mf_sq), 1]])
    B = np.array([(sum_ef / sum_ef_sq), (sum_mf / sum_mf_sq)])
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    return x

# Calculates the average percentage error between A and B.
def average_error_model_a(A, B, x):
    error = 0
    for i, a in enumerate(A):
        a = a[0]
        b = B[i][0]
        if b == 0:
            continue
        error += abs(x*a - b) / b
    error *= 100
    error /= A.shape[0]
    return error

def average_error_model_b(A, M, B, x):
    error = 0
    for i, a in enumerate(A):
        a = a[0]
        mv = M[i]
        b = B[i][0]
        if b == 0:
            continue
        estimate = x[0]*a
#        estimate += x[1]*mv
        error += abs(estimate - b) / b
    error *= 100
    error /= A.shape[0]
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
        print("update type:", update[0][0])
        x = minimize_percentage_error_model_b(A, mv, B)
        print("coefficients:", x)
        if x[1] == 0:
            trained_error = average_error_model_a(A, B, x)
            baseline_error = average_error_model_a(A, B, 1)
        else:
            trained_error = average_error_model_b(A, mv, B, x)
            baseline_error = average_error_model_b(A, mv, B, [1, 1])
        print("trained error:", trained_error)
        print("baseline error:", baseline_error)
        print()

if __name__ == "__main__":
    print_solutions("data2/lowres_17f_target400_data.txt")
