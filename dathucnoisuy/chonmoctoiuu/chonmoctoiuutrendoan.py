import math
import os

# Lấy đường dẫn thư mục chứa chính file script này
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
INP_PATH = os.path.join(CURRENT_DIR, "inp.txt")
OUT_PATH = os.path.join(CURRENT_DIR, "out.txt")

# 1. Đọc số lượng mốc n từ file inp.txt cùng thư mục
try:
    with open(INP_PATH, "r", encoding="utf-8") as f:
        content = f.read().strip()
        n = int(float(content))
except FileNotFoundError:
    print(f"Không tìm thấy file {INP_PATH}! Đặt mặc định n = 5.")
    n = 5
except ValueError:
    print("Dữ liệu trong file inp.txt không hợp lệ! Đặt mặc định n = 5.")
    n = 5

print(f"Số mốc nội suy đọc được từ inp.txt: n = {n}")

# 2. Nhập khoảng [a, b] từ bàn phím
while True:
    try:
        a = float(input("Nhập giá trị a: "))
        b = float(input("Nhập giá trị b: "))
        if a > b:
            # Tự động hoán đổi nếu người dùng nhập a > b
            a, b = b, a
            print(f"Lưu ý: a > b nên đoạn được chuẩn hóa thành [{a}, {b}]")
        break
    except ValueError:
        print("Vui lòng nhập số thực hợp lệ!")

# 3. Tính toán n mốc nội suy Chebyshev tối ưu trên đoạn [a, b]
x_nodes = []
for k in range(1, n + 1):
    # Nghiệm Chebyshev trên đoạn [-1, 1]
    t_k = math.cos(((2 * k - 1) / (2 * n)) * math.pi)
    
    # Chuyển đổi sang đoạn [a, b]
    x_k = ((a + b) / 2) + ((b - a) / 2) * t_k
    x_nodes.append(x_k)

# 4. Sắp xếp các mốc theo thứ tự tăng dần
x_nodes.sort()

# 5. In kết quả ra màn hình
print(f"\n{n} mốc nội suy tối ưu theo thứ tự tăng dần trên [{a}, {b}]:")
for idx, val in enumerate(x_nodes, start=1):
    print(f"x_{idx} = {val:.6f}")

# 6. Ghi kết quả ra file out.txt cùng thư mục
with open(OUT_PATH, "w", encoding="utf-8") as f:
    result = " ".join(f"{val:.6f}" for val in x_nodes)
    f.write(result + "\n")

print(f"\nKết quả đã được ghi vào file out.txt ({OUT_PATH}) thành công!")