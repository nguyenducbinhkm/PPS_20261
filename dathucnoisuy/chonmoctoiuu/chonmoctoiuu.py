import math

# Mở file inp.txt để đọc giá trị n (Lưu ý: trong C++ bạn dùng "w" là chế độ ghi, 
# nhưng để đọc dữ liệu đầu vào thì chuẩn phải là "r")
try:
    with open("inp.txt", "r") as f:
        n = float(f.read().strip())
except FileNotFoundError:
    n = 0.0

n_int = int(n)
x = [0.0] * (n_int + 1)
pi = math.pi

# Tính toán giá trị mảng x
for i in range(n_int + 1):
    x[i] = math.cos(((2 * i + 1) / (2 * n + 2)) * pi) / (2 ** n)

# Ghi kết quả ra file out.txt (tương đương freopen("out.txt", "w", stdout))
with open("out.txt", "w") as f:
    result = " ".join(map(str, x)) + " "
    f.write(result) 