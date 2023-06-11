
# An Example of Solving Bar Element 
# Example is taken from the book: A First Cource in the Finite Element Method, Example 3.1

# E1=207 GP
# A1=654 mm2
# E2=103 GPa
# A2=1290 mm3
# F= 13,345 N

def kael(a,e,l):
    return (a*e)/l

def gaussian_elimination(A, b):
    n = len(A)
    M = A

    i = 0
    for x in M:
        x.append(b[i])
        i += 1

    for k in range(n):
        for i in range(k, n):
            if abs(M[i][k]) > abs(M[k][k]):
                M[k], M[i] = M[i], M[k]
            else:
                pass

        for j in range(k+1, n):
            q = float(M[j][k]) / M[k][k]
            for m in range(k, n+1):
                M[j][m] -= q * M[k][m]

    x = [0 for _ in range(n)]

    x[n-1] = float(M[n-1][n]) / M[n-1][n-1]
    for i in range(n-1, -1, -1):
        z = 0
        for j in range(i+1, n):
            z += float(M[i][j]) * x[j]
        x[i] = float(M[i][n] - z) / M[i][i]
    return x

#Element1
#Element2
matrix=[[1,-1],[-1,1]]
k1 = kael(645*1e-6,207*1e9,0.76)
k2=k1
k3= kael(1290*1e-6,103*1e9,0.76)

kmatrix1 = [[0,0],[0,0]]
for i in range(len(matrix)):
    for j in range(len(matrix)):
        kmatrix1[i][j] = k1*matrix[i][j]

kmatrix2 = [[0,0],[0,0]]
for i in range(len(matrix)):
    for j in range(len(matrix)):
        kmatrix2[i][j] = k2*matrix[i][j]

kmatrix3 = [[0,0],[0,0]]
for i in range(len(matrix)):
    for j in range(len(matrix)):
        kmatrix3[i][j] = k3*matrix[i][j]

K_global = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
for i in range(2):
    for j in range(2):
        K_global[i][j] += kmatrix2[i][j]

for i in range(2):
    for j in range(2):
        K_global[i+1][j+1] += kmatrix3[i][j]

#F=Ku
#Boundary condition
#u1=0, u4=0
#Eliminate the first and fourth columns

K_boundary_rows = [row for i, row in enumerate(K_global) if i not in {0, 3}]

K_boundary = [[col for j, col in enumerate(row) if j not in {0, 3}] for row in K_boundary_rows]

F=[13,345,0]

U= gaussian_elimination(K_boundary,F)

