from math import *
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from sympy import *

np.set_printoptions(
    suppress=True, linewidth=np.inf, precision=5
)  # chỉnh số chữ số sau dấu phẩy


def Input():
  x, y = [], []
  # Tự động lấy đường dẫn tuyệt đối theo thư mục chứa file .py hiện tại
  current_dir = os.path.dirname(os.path.abspath(__file__))
  file_path = os.path.join(current_dir, "input.txt")

  with open(file_path, "r", encoding="utf-8") as f:  # đọc file input
    for line in f.readlines():  # duyệt từng hàng trong file
      parts = line.split()
      if parts:
        val_x = float(parts[0])
        if val_x not in x:
          x.append(val_x)
          y.append(float(parts[1]))
  return x, y


# kết nạp mốc x[i] vào Bảng tỷ hiệu
def BuildBTH(BTH, x, y, i):
  TH_row_i = np.zeros(len(x))
  TH_row_i[0] = y[i]
  for j in range(1, i + 1):
    TH_row_i[j] = (TH_row_i[j - 1] - BTH[j - 1]) / (x[i] - x[i - j])
  return TH_row_i


# kết nạp mốc x[i] vào Bảng tính tích
def BuildBTT(BTT, x, i):
  TT_row_i = np.zeros(len(x))
  TT_row_i[i + 1] = 1  # sử dụng lược đồ Hoocne
  for j in range(i, 0, -1):
    TT_row_i[j] = BTT[j - 1] - x[i] * BTT[j]
  TT_row_i[0] = -x[i] * BTT[0]
  return TT_row_i


# kết nạp mốc x[i] vào đa thức
def fx(BTH, BTT, i, f):
  for j in range(0, i + 1):
    f[j] += BTH[i] * BTT[j]
  return f


# Xây dựng đa thức nội suy Newton tiến
def Newton_interpolation_Forward(BTH, BTT, f, x, y, n):
  for i in range(1, n + 1):
    BTH = BuildBTH(BTH, x, y, i)
    BTT = BuildBTT(BTT, x, i - 1)
    f = fx(BTH, BTT, i, f)
  return f, BTH, BTT


# ước tính f(x0)
def fx0(f, n, x0):
  value = f[n]
  for i in range(n - 1, -1, -1):  # sử dụng lược đồ Hoocne
    value = value * x0 + f[i]
  return value


# Xây dựng đa thức nội suy Newton lùi
def Newton_interpolation_Backward(BTH, BTT, f, x, y, n):
  x.reverse()
  y.reverse()
  BTH[0], f[0] = y[0], y[0]
  return Newton_interpolation_Forward(BTH, BTT, f, x, y, n)


# Thêm mốc nội suy
def add_Newton(BTH, BTT, f, x, y, n):
  BTH = BuildBTH(BTH, x, y, n)
  BTT = BuildBTT(BTT, x, n - 1)
  f = np.hstack((f, np.zeros(1)))  # tạo thêm 1 cột cho f
  f = fx(BTH, BTT, n, f)
  return f, BTH, BTT


# vẽ hình
def plot(x, y, f, n):
  plt.figure()
  plt.title("Newton interpolation")
  plt.scatter(x, y)
  xx = np.linspace(min(x), max(x), 100)
  plt.plot(xx, fx0(f, n, xx), color="blue", label="f_Newton")
  plt.legend()
  current_dir = os.path.dirname(os.path.abspath(__file__))
  plt.savefig(os.path.join(current_dir, "newton.png"))
  plt.close()


def main():
  x, y = Input()
  n = len(x) - 1
  x0 = float(input("Mời nhập giá trị cần tính: "))

  print("Các mốc nội suy: ")
  print(x)

  # khởi tạo ma trận không 1 chiều
  BTH, BTT, f = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
  BTH[0], BTT[0], f[0] = y[0], 1, y[0]

  f, BTH, BTT = Newton_interpolation_Forward(BTH, BTT, f, x, y, n)

  value = fx0(f, n, x0)
  print("Hệ số đa thức nội suy Newton là: (bắt đầu từ hệ số tự do)")
  print(f)
  print("Giá trị cần tính tại", x0, " là:", value)

  plot(x, y, f, n)

  try:
    k_input = input("Nhập số lượng thêm mốc nội suy (nhấn Enter để bỏ qua): ")
    if k_input.strip() == "":
      return
    k = int(k_input)
    if k <= 0:
      return
    for i in range(k):
      xt, yt = input(f"Nhập mốc và giá trị nội suy thứ {i+1} (cách nhau khoảng trắng): ").split()
      if float(xt) not in x:
        n += 1
        x.append(float(xt))
        y.append(float(yt))
        f, BTH, BTT = add_Newton(BTH, BTT, f, x, y, n)

    value = fx0(f, n, x0)
    print("Hệ số đa thức nội suy Newton mới là: (bắt đầu từ hệ số tự do)")
    print(f)
    print("Giá trị cần tính tại", x0, " là:", value)
    plot(x, y, f, n)
  except Exception as e:
    print("Kết thúc chương trình hoặc bỏ qua phần thêm mốc.")


if __name__ == "__main__":
  main()