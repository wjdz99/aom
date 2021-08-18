import numpy as np

# Uses least squares regression to find the solution
# when there is one unknown variable.
def print_lstsq_solution(A, B):
    A_inv = np.linalg.pinv(A)
    x = np.matmul(A_inv, B)
    return x

# Uses the pseudoinverse matrix to find the solution
# when there are two unknown variables.
def print_pinv_solution(A, mv, B):
    new_A = np.concatenate((A, mv), axis=1)
    print(new_A.shape)
    new_A_inv = np.linalg.pinv(new_A)
    new_x = np.matmul(new_A_inv, B)
    return new_x

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

# Calculates the average percentage error between A and B.
def average_error(A, B, x):
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
        x = minimize_percentage_error_model_a(A, B)
        print("coefficient:", x)
        trained_error = average_error(A, B, x)
        baseline_error = average_error(A, B, 1)
        print("trained error:", trained_error)
        print("baseline error:", baseline_error)
        print()

if __name__ == "__main__":
    print_solutions("data2/lowres_64f_target150_data.txt")
