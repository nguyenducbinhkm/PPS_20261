import sys
import os

# Khai báo cấu trúc dữ liệu tương đương struct Data trong C++
class Data:

    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y


# Biến toàn cục
n = 0
dw = [0.0] * 10012
arr_wi = []  # Hệ số của (x-x_1)(x-x_2)...(x-x_n)
arr_Lagrange = []  # Hệ số của đa thức ns Lagrange
a = [Data() for _ in range(10012)]


def Print(v):
    print(f"In cac he so tu x^0 den x^{len(v) - 1}")
    print(" ".join(map(str, v)))
    print()


def hoocnerNhan(arr, hs, xk):  # arr nhân (hs * x + xk)
    ret_size = len(arr) + 1
    ret = [0.0] * ret_size

    if len(arr) == 0:
        return [xk, hs]

    for i in range(ret_size - 1, -1, -1):
        tmp = 0.0
        if i > 0:
            tmp += arr[i - 1] * hs
        if i < len(arr):
            tmp += arr[i] * xk
        ret[i] = tmp

    return ret


def hoocnerChia(arr, hs, xk):  # arr chia (hs * x + xk)
    ret_size = len(arr) - 1
    ret = [0.0] * ret_size

    nho = 0.0
    for i in range(len(arr) - 1, 0, -1):
        val = (arr[i] - nho) / hs  # hệ số x^(i-1)
        nho = val * xk
        ret[i - 1] = val

    return ret


def Get_Lagrange(arr_wi, pos):
    global arr_Lagrange
    ret = hoocnerChia(arr_wi, 1.0, -a[pos].x)

    val = a[pos].y / dw[pos]
    for i in range(len(arr_Lagrange)):
        arr_Lagrange[i] += ret[i] * val


def Init():
    global arr_wi, arr_Lagrange, n

    # Sắp xếp theo các mốc nội suy tăng dần
    a[: n + 1] = sorted(a[: n + 1], key=lambda item: item.x)

    print("Sap xep theo cac moc noi suy tang dan")
    for i in range(n + 1):
        print(f"{a[i].x} {a[i].y}")
    print()

    # Tính dw[i] = w'(x_i)
    for i in range(n + 1):
        dw[i] = 1.0
        for j in range(n + 1):
            if i != j:
                dw[i] *= a[i].x - a[j].x

    # arr_wi là hệ số của (x-x_1)(x-x_2)...(x-x_n)
    for i in range(n + 1):
        arr_wi = hoocnerNhan(arr_wi, 1.0, -a[i].x)

    # arr_Lagrange chứa hệ số của đa thức nội suy Lagrange
    arr_Lagrange = [0.0] * (len(arr_wi) - 1)
    for i in range(n + 1):
        Get_Lagrange(arr_wi, i)

    Print(arr_Lagrange)


def Get_Value_y(value_x):
    value_y = 0.0
    powx = 1.0
    for i in range(len(arr_Lagrange)):
        value_y += arr_Lagrange[i] * powx
        powx *= value_x
    return value_y


def main():
    global n

    # Tự động lấy thư mục chứa file nslagrange.py hiện tại
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Ghép nối chính xác đường dẫn file inp và out cùng thư mục với file py
    input_path = os.path.join(current_dir, "nsLagrange.inp")
    output_path = os.path.join(current_dir, "nsLagrange.out")

    try:
        sys.stdin = open(input_path, "r", encoding="utf-8")
        sys.stdout = open(output_path, "w", encoding="utf-8")
    except FileNotFoundError:
        print(f"Không tìm thấy file tại: {input_path}")
        return

    input_data = sys.stdin.read().split()
    if not input_data:
        return

    iterator = iter(input_data)

    n = int(next(iterator))  # số lượng mốc nội suy ban đầu
    n -= 1

    for i in range(n + 1):
        a[i].x = float(next(iterator))
        a[i].y = float(next(iterator))

    Init()

    nq = int(next(iterator))  # số lượng mốc nội suy cần tính
    for _ in range(nq):
        value_x = float(next(iterator))
        print(f"Gia tri da thuc ns Lagrange tai x = {value_x} bang ", end="")
        print(f"{Get_Value_y(value_x):.10f}")

if __name__ == "__main__":
    main()