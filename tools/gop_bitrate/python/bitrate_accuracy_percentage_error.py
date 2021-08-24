import numpy as np

def lstsq_solution(A, B):
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    return x[0][0]

# Finds the coefficient to multiply A by to minimize
# the percentage error between A and B.
def minimize_percentage_error_model_a(A, B):
    R = np.divide(A, B)
    num = 0
    den = 0
    best_x = 0
    best_error = 100
    for r_i in R:
        num += r_i
        den += r_i**2
        curr_x = (num/den)[0]
#        print("test", curr_x)
        curr_error = average_error_model_a(A, B, curr_x)
#        print()

        if curr_error < best_error:
            best_error = curr_error
            best_x = curr_x
    if den == 0:
        x = 0
    else:
        x = (num / den)[0]
#    print(x)
    return best_x

def minimize_percentage_error_model_b(r_e, r_m, r_f):
    r_ef = np.divide(r_e, r_f)
    r_mf = np.divide(r_m, r_f)
    sum_ef = np.sum(r_ef)
    sum_ef_sq = np.sum(np.square(r_ef))
    sum_mf = np.sum(r_mf)
    sum_mf_sq = np.sum(np.square(r_mf))
    sum_ef_mf = np.sum(np.multiply(r_ef, r_mf))

    # Divides x by y. If y is zero, returns 0.
    divide = lambda x, y : 0 if y == 0 else x / y

    # Set up and solve the matrix equation
    A = np.array([[1, divide(sum_ef_mf, sum_ef_sq)],[divide(sum_ef_mf, sum_mf_sq), 1]])
    B = np.array([divide(sum_ef, sum_ef_sq), divide(sum_mf, sum_mf_sq)])
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    return x

# Calculates the average percentage error between A and B.
def average_error_model_a(A, B, x):
    actual_error = 0
    estimated_error = 0
    n = 0
    for i, a in enumerate(A):
        a = a[0]
        b = B[i][0]
        if b == 0:
            continue
        actual_error_i = (abs(a-b)/b)*100
        estimated_error_i = (abs(x*a-b)/b)*100
        n += 1
        actual_error += actual_error_i
        estimated_error += estimated_error_i
#        print(actual_error/n, estimated_error/n, "\t", actual_error_i, estimated_error_i, "\t", a, x*a, b)
    actual_error /= n
    estimated_error /= n
#    print(actual_error, estimated_error)
    return estimated_error

def average_error_model_b(A, M, B, x):
    error = 0
    for i, a in enumerate(A):
        a = a[0]
        mv = M[i]
        b = B[i][0]
        if b == 0:
            continue
        estimate = x[0]*a
        estimate += x[1]*mv
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
#        if update[0][0] != 3:
#            continue

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
        xa = minimize_percentage_error_model_a(A, B)
        xb = minimize_percentage_error_model_b(A, mv, B)
        print("Model A coefficients:", xa, " | Model B coefficients:", xb)
        error_a = average_error_model_a(A, B, xa)
        error_b = average_error_model_b(A, mv, B, xb)
        baseline_error_a = average_error_model_a(A, B, 1)
        baseline_error_b = average_error_model_b(A, mv, B, [1, 1])
        print("error a:", error_a, " | error b:", error_b)
        print("baseline error a:", baseline_error_a, "baseline error b:", baseline_error_b)
        print()

if __name__ == "__main__":
    print_solutions("data2/lowres_17f_target400_data.txt")
