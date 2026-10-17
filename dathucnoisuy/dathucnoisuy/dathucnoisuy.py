import sys

# Khai báo các biến toàn cục tương đương C++
n_val = -1
N = 0
x = [0.0] * 10001
y = [0.0] * 10001
X = [0.0] * 10001
Y = [0.0] * 10001
kq = [0.0] * 10001
Van = [[0.0] * 1001 for _ in range(1001)]


def kiemtratontai():
    global n_val
    for i in range(N + 1):
        check = True
        for j in range(i + 1, N + 1):
            if X[i] == X[j]:
                check = False
        if check:
            n_val += 1
            x[n_val] = X[i]
            y[n_val] = Y[i]


def taoVandemonde():
    for i in range(n_val + 1):
        Van[i][0] = 1.0
        for j in range(1, n_val + 1):
            tg = x[i]
            Van[i][j] = 1.0
            for k in range(1, j + 1):
                Van[i][j] *= tg


def giaimatranVan():
    for j in range(n_val):
        for i in range(j + 1, n_val + 1):
            temp = Van[i][j] / Van[j][j]
            for k in range(n_val + 2):
                Van[i][k] -= Van[j][k] * temp

    for i in range(n_val, -1, -1):
        s = 0.0
        for j in range(n_val, i - 1, -1):
            s += Van[i][j] * kq[j]
        kq[i] = (Van[i][n_val + 1] - s) / Van[i][i]


def tinhgiatri(c):
    b = [0.0] * 10001
    b[n_val] = kq[n_val]
    for i in range(n_val - 1, -1, -1):
        b[i] = b[i + 1] * c + kq[i]
    gt = b[0]
    return gt


def main():
    global N
    # Đọc từ file inp.txt và ghi ra out.txt (tương đương freopen)
    try:
        sys.stdin = open("inp.txt", "r")
        sys.stdout = open("out.txt", "w")
    except FileNotFoundError:
        pass

    # Đọc N và các cặp tọa độ (X[i], Y[i])
    input_data = sys.stdin.read().split()
    if not input_data:
        return

    iterator = iter(input_data)
    N = int(next(iterator))

    for i in range(N + 1):
        X[i] = float(next(iterator))
        Y[i] = float(next(iterator))

    kiemtratontai()

    print("-" * 40)
    print("Cac moc noi suy ban dau la:")
    for i in range(n_val + 1):
        print(f"{x[i]:8.7g}    {y[i]:7.7f}  ")

    print("-" * 40)
    taoVandemonde()
    print("Ma tran Vandemonde tao boi cac moc noi suy:")
    for i in range(n_val + 1):
        Van[i][n_val + 1] = y[i]

    for i in range(n_val + 1):
        row_str = "".join(f"{Van[i][j]:12.7g}" for j in range(n_val + 1))
        print(row_str)

    print("-" * 40)
    giaimatranVan()
    print("Ma tran mo rong cua hpt duoc chuyen thanh ma tran bac thang nhu sau:")
    for i in range(n_val + 1):
        row_str = "".join(f"{Van[i][j]:8.7g}  " for j in range(n_val + 1))
        print(row_str)

    print("-" * 40)
    print("Da thuc noi suy la: f(x) = ", end="")
    for i in range(n_val + 1):
        if i == 0:
            print(f"a{i} + ", end="")
        elif i == 1:
            print(f"a{i}x + ", end="")
        elif i == n_val:
            print(f"a{i}x^{i}", end="")
        else:
            print(f"a{i}x^{i} + ", end="")
    print()

    for i in range(n_val + 1):
        print(f"x{i}={kq[i]:8.8g}  ", end="")
    print()

    print("-" * 40)
    q = float(next(iterator))
    print(f"Tinh gia tri cua bieu thuc tai diem x = {q} la:")
    print(f"{tinhgiatri(q):.10f}")


if __name__ == "__main__":
    main()