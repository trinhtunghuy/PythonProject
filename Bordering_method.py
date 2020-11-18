import time
import numpy as np
np.set_printoptions(suppress=True, linewidth=np.inf, precision=12)
print("========================================================================")
print("Chương Trình Tìm nghịch đảo của ma trận bằng phương pháp viền quanh")
print("========================================================================")
a = np.loadtxt("test.txt",dtype='float', delimiter=' ')
print(a)
det = np.linalg.det(a)
if det == 0:
    print("Ma trận không khả nghịch!!!")
    quit()
n = len(a)
b = np.transpose(a)
t = b.dot(a)
print("Ma trận vừa nhập là:")
print(a)
print("Nhân A^T với A ta được:")
print(t)
print("========================================================================")
def bordering(a, n):# Nhap ma tran a va kich thuoc cua ma tran

    if n == 2: #Trường hợp chặn đệ quy
        inv = np.zeros([2,2]) #Khởi tạo ma trận nghịch đảo alpha^-1
        inv[0,0] = a[1,1]/(a[0,0] * a[1,1] - a[1,0]*a[0, 1])
        inv[0,1] = -a[0 ,1]/(a[0,0] * a[1,1] - a[1,0]*a[0, 1])
        inv[1,0] = -a[1,0]/(a[0,0] * a[1,1] - a[1,0]*a[0, 1])
        inv[1,1] = a[0,0]/(a[0,0] * a[1,1] - a[1,0]*a[0, 1])
        print("A_11^-1=")
        print(inv)
        return inv
    a_11 = np.zeros([n - 1, n - 1])  # Khởi tạo các ma trận con
    a_12 = np.zeros([n - 1, 1])
    a_21 = np.zeros([1, n - 1])
    a_22 = a[n - 1, n - 1]
    for i in range(0, n - 1):
        for j in range(0, n - 1):
            a_11[i, j] = a[i, j]
    print("A_11:")
    print(a_11)
    for i in range(0, n - 1):
        a_12[i, 0] = a[i, n - 1]
        a_21[0, i] = a[n - 1, i]
    print("A_12:")
    print(a_12)
    print("A_21:")
    print(a_21)
    print("A_22:")
    print(a_22)
    print("====================================================")
    if n>2: #Nếu chưa gặp trường hợp chặn
        a_111 = bordering(a_11, n-1)
        X = a_111.dot(a_12)
        Y = a_21.dot(a_111)
        theta = a_22 - Y.dot(a_12)
        print("X=")
        print(X)
        print("Y=")
        print(Y)
        print("Theta = ")
        print(theta)
        x = X.dot(Y)
        for i in range(0, n - 1):
            for j in range(0, n - 1):
                a[i, j] = a_111[i, j] + x[i, j] / theta
                a[n - 1, i] = Y[0, i] / (-theta)
                a[i, n - 1] = -X[i, 0] / theta
        a[n - 1, n - 1] = 1 / theta
        print("A_11^-1")
        print(a)
    return a
print("========================================================================")
st = time.time()
s = bordering(t, n)
et = time.time()
print("========================================================================")
print("Ma trận nghịch đảo là:")
print(s.dot(b))
print("========================================================================")
print((et-st)*1000,'ms')
print("========================================================================")
print("Kiểm tra lại kết quả:")
print(a.dot(s.dot(b)))
