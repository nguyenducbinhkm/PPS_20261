import sys

# Khai báo các biến toàn cục
a = [0.0] * 10005
x = 0.0
n = 0
gt = 0.0


# Lớp giả lập mảng C++ (hỗ trợ cả chỉ số âm như b[-1])
class CppArray:

    def __init__(self):
        self.data = {}

    def __getitem__(self, index):
        return self.data.get(index, 0.0)

    def __setitem__(self, index, value):
        self.data[index] = value


b = CppArray()

# Tinh gia tri cua da thuc tai 1 diem
def tinhgiatri(c):
    global gt
    b[n] = a[n]
    for i in range(n - 1, -1, -1):
        b[i] = b[i + 1] * c + a[i]
    gt = b[0]
    return gt


# Tinh gia tri dao ham cua da thuc tai 1 diem
def tinhdaoham(c):
    global gt
    d = {}
    for i in range(n, -1, -1):
        d[i - 1] = i * a[i]
    b[n - 1] = d[n - 1]
    for i in range(n - 2, -1, -1):
        b[i] = b[i + 1] * c + d[i]
    gt = b[0]
    return gt


# Tim da thuc thuong va phan du khi lay da thuc chia cho don thuc
def chiadathuc(c):
    b[n] = a[n]
    for i in range(n - 1, -1, -1):
        b[i] = b[i + 1] * c + a[i]
    out = []
    for i in range(n, 0, -1):
        out.append(str(b[i]))
    out.append(str(b[0]))
    print(" ".join(out), end="")


# Tim da thuc tich 
def nhandathuc(c):
    b[n + 1] = a[n]
    b[n] = a[n]
    b[-1] = -a[0] * c
    for i in range(n - 1, -1, -1):
        b[i] = a[i] - a[i + 1] * c
    out = []
    for i in range(n, -2, -1):
        out.append(str(b[i]))
    print(" ".join(out), end="")


def main():
    global n, x, a

    # Mở file inp.txt để đọc và out.txt để ghi kết quả
    try:
        sys.stdin = open("inp.txt", "r")
        sys.stdout = open("out.txt", "w")
    except FileNotFoundError:
        pass

    input_data = sys.stdin.read().split()
    if not input_data:
        return

    iterator = iter(input_data)
    n = int(next(iterator))

    for i in range(n, -1, -1):
        a[i] = float(next(iterator))

    x = float(next(iterator))
    t = int(next(iterator))

    if t == 1:
        print(f"Gia tri cua da thuc f(x) tai diem x = {x} la {tinhgiatri(x)}")
    elif t == 2:
        print(
            f"Gia tri cua dao ham cua da thuc tai diem x = {x} la {tinhdaoham(x)}"
        )
    elif t == 3:
        chiadathuc(x)
    elif t == 4:
        nhandathuc(x)


if __name__ == "__main__":
    main()