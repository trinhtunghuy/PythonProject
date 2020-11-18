from cmath import sqrt
import numpy as np
def cholesky():
    S=[[0]*n for i in range(n)]
    for i in range(n):
        for j in range(i + 1):
            tổng = sum(S[k][i] * S[k][j] for k in range(j))
            if abs(tổng - A1[i][i]) == 0: break
            if i == j: S[j][i] = sqrt(A1[i][i] - tổng)
            else: S[j][i] = (A1[j][i] - tổng)/S[j][j]
        if S[i][i] == 0: break
    return S
def cholesky3(): #Hàm giải phương trình BY = E
    Y = np.zeros([n, n],dtype="complex_")
    for i in range(n-1, -1, -1):
        for j in range(i,n):
            if j == i: Y[i,j] = 1/S[i,j]
            if j>i:
                Y[i, j] = (B1[i,j]-sum(S[i][k]*Y[k,j] for k in range(i+1,n)))/S[i][i]
    return Y
#Nhập ma trận từ file
A = np.loadtxt("test.txt",dtype='complex', delimiter=' ')
k=0
print('The input matrix is:')
print(A)
n = len(A)
B = np.zeros([n,n])
for i in range(n):
    B[i,i] = 1.0
print('')
print('Tìm ma trận A1 đối xứng:')
m = 1
for i in range(n):
    for j in range(n):
        if A[i][j] == A[j][i]: continue
        else: m *= 0
if m == 1:
        print('A là ma trận vuông đối xứng có thể khai triển theo Cholesky')
        A1=A
        B1=B
else:
        print('A không đối xứng')
        print('Để giải được theo Cholesky ta cần nhân cả 2 vế với A^t:')
        print("A*A^t=")
        A1 = A.dot(np.transpose(A))
        B1 = B
        k = 1
print(A1)
S = np.array(cholesky())
if S[n-1][n-1] == 0:
    print('Ma trận A1 có định thức bằng 0 nên không thể khai triển Cholesky')
    quit()
else:
    print('')
    print('Phân tích A1 theo Cholesky: A1=S^t*S')
    print('P=')
    print(str(np.array(S)))
    print('Q=')
    s = np.transpose(S)
    print(s)
Y = cholesky3()
print("==============================================================")
print("P^-1=")
print(Y)
Y = Y.dot(np.transpose(Y))
print()
print("=============================================================")
print('The inverse matrix is:')
if k == 0:
    print(Y)
if k == 1:
    Y = np.transpose(A).dot(Y)
    print(Y)
print('===============================================================')
print("Test if the solution is true:")
print("A*A^-1=")
print(A.dot(Y))