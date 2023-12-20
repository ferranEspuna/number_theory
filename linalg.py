import numpy as np


def find_first_nonzero(row):
    for i, val in enumerate(row):
        if val != 0:
            return i, val
    return None


# Finds the row-reduced echelon form of a matrix over a field
# The matrix is modified in-place
# Returns the rank of the matrix
def reduce(matrix, char):
    n, m = matrix.shape
    for i in range(n):

        first_own_nonzero = find_first_nonzero(matrix[i])
        best_other_nonzero = None
        for j in range(i + 1, m):
            other_nonzero = find_first_nonzero(matrix[j])
            if other_nonzero is None:
                continue
            if best_other_nonzero is None or other_nonzero[0] < best_other_nonzero[1]:
                best_other_nonzero = (j, *other_nonzero)

        if best_other_nonzero is None:
            if first_own_nonzero is None:
                return i
            else:
                pos, val = first_own_nonzero
                val = int(val)
                inv = pow(val, -1, char)
                matrix[i] = (matrix[i] * inv) % char
                return i + 1

        j, k, val = best_other_nonzero
        if first_own_nonzero is None or k < first_own_nonzero[0]:
            matrix[[i, j]] = matrix[[j, i]]
            first_own_nonzero = best_other_nonzero[1:]

        pos, val = first_own_nonzero
        val = int(val)
        inv = pow(val, -1, char)
        matrix[i] = (matrix[i] * inv) % char

        for j in range(i + 1, m):
            other_nonzero = find_first_nonzero(matrix[j])
            if other_nonzero is None:
                continue
            if other_nonzero[0] == first_own_nonzero[0]:
                matrix[j] = (matrix[j] - matrix[i] * other_nonzero[1]) % char

    return n


# Finds a non-trivial (different to 0, and to (k, 0,..., 0) solution to the equation Ax = 0
# where A is a matrix over a field that is in reduced row echelon form
def find_some_non_trivial_in_reduced(matrix, char):
    n, m = matrix.shape
    i = n - 1
    j = n
    candidate = np.zeros(m, dtype=int)
    found = False
    while i >= 0:

        t = find_first_nonzero(matrix[i])

        if i == 0 and not found and t is None:
            t = (0, 1)

        if t:

            new_j, val = t
            jump = j - new_j
            assert jump > 0, 'Matrix is not in reduced form!'
            assert val == 1, 'Matrix is not in reduced form!'
            j = new_j

            if found:
                defect = np.dot(matrix[i], candidate) % char
                if defect != 0:
                    candidate[j] = char - defect

            else:

                if jump > 1:
                    candidate[j + 1] = 1
                    candidate[j] = (-matrix[i, j + 1]) % char
                    found = True

        i -= 1

    if not found:
        return None

    return candidate


def find_faster(matrix, char):

    i = 0
    while matrix[i, i+1] != 0:
        i += 1

    lvec = i+2
    sol = np.zeros(lvec, dtype=int)
    sol[i+1] = 1

    while i > 0:
        i -= 1

        res = np.dot(matrix[i, :lvec], sol) % char
        sol[i+1] = (-res) % char

    return sol