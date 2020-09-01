# Implements the simplex algorithm using standard python3 libraries

# if |a-b| < tol then a and b are treated as equal
tol = 1e-9

class Set():
    # data structure for tracking the basic/nonbasic variables on each iteration of simplex
    def __init__(self, seq):
        self.idx = {var:idx for idx, var in enumerate(seq)}
        self.var = list(seq)

    def replace(self, var_1, var_2):
        assert var_1 in self.idx
        i = self.idx.pop(var_1)
        assert self.var[i] == var_1
        self.var[i] = var_2
        self.idx[var_2] = i
    
    def add(self, new):
        i = len(self.var)
        self.var.append(new)
        assert self.var[i] == new
        self.idx[new] = i
        

def simplex(A, b, c):
    assert len(A) == len(b) and all([len(row) == len(c) for row in A]), "incorrect dimensions for input"
    if all(x < tol for x in b) and all(x < tol for x in c):
        return 0, [0]*len(c)
    L_aux = initialize_simplex(A, b, c)
    if L_aux:
        N, B, A, b, c, v = L_aux
    else:
        print("No initial feasible solution found. LP has no solutions.")
        return list()
    # check indices in reverse on every 5th iteration to avoid cycling
    counter = 1
    while any([x > tol for x in c]):
        if counter % 5:
            j, e = next((j, N.var[j]) for j in range(len(c)) if c[j] > tol)
        else:
            j, e = next((j, N.var[j]) for j in reversed(range(len(c))) if c[j] > tol)
        counter += 1
        temp = [b[i]/A[i][j] if A[i][j] > tol else float("inf") for i in range(len(b))]
        i, val = min(enumerate(temp), key = lambda x: x[1])
        if val == float("inf"):
            print("No bounded optimal solution.")
            return []
        else:
            l = B.var[i]
            N, B, A, b, c, v = pivot(N, B, A, b, c, v, l, e)
    print("Bounded solution.")
    return [(b[B.idx[i]] if (i in B.var) else 0) for i in range(len(c))]

def pivot(N, B, A, b, c, v, l, e):
    # performs the pivot operation used in each iteration of simplex
    row, col = B.idx[l], N.idx[e]
    m, n = len(b), len(c)
    div = A[row][col]
    b[row] /= div
    A[row][col] = 1.0
    A[row] = [x/div for x in A[row]]
    # compute coefficients not in pivot row/col
    for i in range(m):
        if i == row: continue
        b[i] -= A[i][col]*b[row]
        for j in range(n):
            if j == col: continue
            A[i][j] -= A[i][col]*A[row][j]
        A[i][col] /= -div
    v += c[col]*b[row]
    for j in range(n):
        if j == col: continue
        c[j] -= c[col]*A[row][j]
    c[col] /= -div
    N.replace(e,l)
    B.replace(l,e)
    return N, B, A, b, c, v

def initialize_simplex(A,b,c):
    # m equations, n variables
    m, n = len(b), len(c)
    assert len(A) == m and all(len(x) == n for x in A), "incorrect dimensions"
    N = Set(range(n))
    B = Set(range(n,n+m))
    k = min(range(len(b)), key = lambda x: b[x])
    if b[k] > tol:
        return N, B, A, b, c, 0
    # form the auxilliary LP
    obj = [x for x in c]
    c = [0]*(n) + [-1]
    A = [row + [-1] for row in A]
    N.add(n+m)
    N, B, A, b, c, v = pivot(N, B, A, b, c, 0, n+k, n+m)
    # use simplex loop to find optimal solution to aux LP
    counter = 1
    while any([x > tol for x in c]):
        if counter % 5:
            j, e = next((j, N.var[j]) for j in range(len(c)) if c[j] > tol)
        else:
            j, e = next((j, N.var[j]) for j in reversed(range(len(c))) if c[j] > tol)
        counter += 1
        temp = [b[i]/A[i][j] if A[i][j] > tol else float("inf") for i in range(len(b))]
        i, val = min(enumerate(temp), key = lambda x: x[1])
        if val == float("inf"):
            return list()
        else:
            l = B.var[i]
            N, B, A, b, c, v = pivot(N, B, A, b, c, v, l, e)
    x = [(b[B.idx[i]] if (i in B.idx) else 0) for i in range(n+m+1)]
    if abs(x[n+m]) < tol:
        if n+m in B.idx:
            # perform one pivot to make auxilliary variable nonbasic
            i, l = B.idx[n+m], n+m
            e = next(e for e, j in N.idx.items() if abs(A[i][j]) > tol) 
            N, B, A, b, c, v = pivot(N, B, A, b, c, v, l, e)
        # remove column corresponding to auxilliary variable
        col = N.idx[n+m]
        A = [row[:col]+row[col+1:] for row in A]
        c = [0]*(n)
        v = 0
        N = Set(N.var[:col]+N.var[col+1:])
        # replace basic variables in objective function with their constraints
        assert len(obj) == len(c)
        for i, val in enumerate(obj):
            if abs(val) < tol: continue
            elif i in B.var:
                index = B.idx[i]
                temp = [-val*x for x in A[index]]
                for j in range(len(c)):
                    c[j] += temp[j]
                v += val*b[index]
            else:
                index = N.idx[i]
                c[index] += val
        return N, B, A, b, c, v
    else:
        # no feasible solution found
        return list()
