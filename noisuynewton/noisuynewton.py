from math import *
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(
    suppress=True, linewidth=np.inf, precision=5
)  # Chỉnh số chữ số sau dấu phẩy


def Input():
    x, y = [], []
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "input.txt")

    if not os.path.exists(file_path):
        print(f"Lỗi: Không tìm thấy file {file_path}")
        return x, y

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f.readlines():
            parts = line.split()
            if parts:
                val_x = float(parts[0])
                if val_x not in x:
                    x.append(val_x)
                    y.append(float(parts[1]))
    return x, y


# Tính bảng tỷ hiệu đầy đủ để lấy các hệ số Newton
def ComputeDividedDifferences(x, y):
    n = len(x)
    table = np.zeros((n, n))
    for i in range(n):
        table[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (
                x[i + j] - x[i]
            )

    return table


# Định dạng đa thức dạng Newton: P(x) = f[x0] + f[x0,x1](x-x0) + ...
def format_newton_form(x, coeffs):
    terms = []
    for i, c in enumerate(coeffs):
        if abs(c) < 1e-9 and i > 0:
            continue
        c_str = f"{c:g}" if abs(c - round(c)) < 1e-5 else f"{c:.5f}".rstrip("0").rstrip(".")

        if i == 0:
            terms.append(f"{c_str}")
        else:
            factors = "".join(
                [
                    f"(x - {x[j]:g})" if x[j] >= 0 else f"(x + {-x[j]:g})"
                    for j in range(i)
                ]
            )
            sign = " + " if c > 0 else " - "
            val = abs(c)
            val_str = f"{val:g}" if abs(val - round(val)) < 1e-5 else f"{val:.5f}".rstrip("0").rstrip(".")
            term = factors if val == 1 else f"{val_str}*{factors}"
            terms.append(f"{sign}{term}")

    return "".join(terms) if terms else "0"


# Định dạng đa thức dạng đại số chuẩn thu gọn: a_n*x^n + ... + a_0
def format_standard_poly(f):
    terms = []
    degree = len(f) - 1
    # f[i] là hệ số của x^i, ta duyệt từ bậc cao nhất về 0
    for p in range(degree, -1, -1):
        coef = f[p]
        if abs(coef) < 1e-9:
            continue
        val = abs(coef)
        sign = " - " if coef < 0 else " + "
        val_str = f"{val:g}" if abs(val - round(val)) < 1e-5 else f"{val:.5f}".rstrip("0").rstrip(".")

        if p == 0:
            term = f"{val_str}"
        elif p == 1:
            term = "x" if val == 1 else f"{val_str}*x"
        else:
            term = f"x^{p}" if val == 1 else f"{val_str}*x^{p}"

        if not terms:
            terms.append(f"-{term}" if coef < 0 else term)
        else:
            terms.append(f"{sign}{term}")

    return "".join(terms) if terms else "0"


# Kết nạp mốc x[i] vào Bảng tỷ hiệu (lũy kế)
def BuildBTH(BTH, x, y, i):
    TH_row_i = np.zeros(len(x))
    TH_row_i[0] = y[i]
    for j in range(1, i + 1):
        TH_row_i[j] = (TH_row_i[j - 1] - BTH[j - 1]) / (x[i] - x[i - j])
    return TH_row_i


# Kết nạp mốc x[i] vào Bảng tính tích
def BuildBTT(BTT, x, i):
    TT_row_i = np.zeros(len(x))
    TT_row_i[i + 1] = 1  # Sử dụng lược đồ Hoocne
    for j in range(i, 0, -1):
        TT_row_i[j] = BTT[j - 1] - x[i] * BTT[j]
    TT_row_i[0] = -x[i] * BTT[0]
    return TT_row_i


# Kết nạp mốc x[i] vào đa thức
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


# Ước tính f(x0) theo lược đồ Horner
def fx0(f, n, x0):
    value = f[n]
    for i in range(n - 1, -1, -1):
        value = value * x0 + f[i]
    return value


# Thêm mốc nội suy
def add_Newton(BTH, BTT, f, x, y, n):
    BTH = BuildBTH(BTH, x, y, n)
    BTT = BuildBTT(BTT, x, n - 1)
    f = np.hstack((f, np.zeros(1)))
    f = fx(BTH, BTT, n, f)
    return f, BTH, BTT


# Xuất kết quả ra file outputnewton.txt
def ExportOutput(x, y, f, x0=None, value=None, is_added=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(current_dir, "outputnewton.txt")
    mode = "a" if is_added else "w"

    dd_table = ComputeDividedDifferences(x, y)
    newton_coeffs = dd_table[0, :]  # Các hệ số tỷ hiệu cấp 0, 1, ..., n

    with open(out_file, mode, encoding="utf-8") as f_out:
        if not is_added:
            f_out.write("=" * 70 + "\n")
            f_out.write("         KẾT QUẢ ĐA THỨC NỘI SUY NEWTON TIẾN\n")
            f_out.write("=" * 70 + "\n\n")
        else:
            f_out.write("\n" + "=" * 70 + "\n")
            f_out.write("     KẾT QUẢ SAU KHI BỔ SUNG THÊM MỐC NỘI SUY\n")
            f_out.write("=" * 70 + "\n\n")

        # 1. Bảng mốc
        f_out.write("1. BẢNG MỐC NỘI SUY:\n")
        f_out.write(f"{'i':<6}{'x_i':<15}{'y_i':<15}\n")
        f_out.write("-" * 36 + "\n")
        for i in range(len(x)):
            f_out.write(f"{i:<6}{x[i]:<15.5f}{y[i]:<15.5f}\n")
        f_out.write("\n")

        # 2. Các hệ số tỷ hiệu Newton
        f_out.write("2. CÁC HỆ SỐ TỶ HIỆU f[x_0, ..., x_k]:\n")
        for k, c in enumerate(newton_coeffs):
            f_out.write(f"   c_{k} = f[x_0...x_{k}] = {c:.6f}\n")
        f_out.write("\n")

        # 3. Biểu diễn đa thức
        f_out.write("3. BIỂU THỨC ĐA THỨC NỘI SUY NEWTON:\n")
        f_out.write(f"   a) Dạng Newton chính tắc:\n")
        f_out.write(f"      P(x) = {format_newton_form(x, newton_coeffs)}\n\n")

        f_out.write(f"   b) Dạng khai triển thu gọn (lũy thừa của x):\n")
        f_out.write(f"      P(x) = {format_standard_poly(f)}\n\n")

        f_out.write(f"   c) Mảng hệ số từ x^0 đến x^{len(f)-1}:\n")
        f_out.write(f"      {list(np.round(f, 6))}\n\n")

        # 4. Giá trị tại x0 (nếu có)
        if x0 is not None and value is not None:
            f_out.write("4. GIÁ TRỊ NỘI SUY TẠI ĐIỂM CẦN TÍNH:\n")
            f_out.write(f"   P({x0}) = {value:.6f}\n")

    print(f"\n>>> Đã xuất kết quả vào file: {out_file}")


# Vẽ hình
def plot(x, y, f, n):
    plt.figure()
    plt.title("Newton interpolation")
    plt.scatter(x, y, color="red", label="Data points")
    xx = np.linspace(min(x), max(x), 200)
    plt.plot(xx, fx0(f, n, xx), color="blue", label="P_Newton(x)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(current_dir, "newton.png")
    plt.savefig(img_path, dpi=300)
    plt.close()
    print(f">>> Đã lưu đồ thị vào file: {img_path}")


def main():
    x, y = Input()
    if not x:
        return

    n = len(x) - 1
    print("Các mốc nội suy đọc được từ input.txt:")
    for i in range(len(x)):
        print(f"  Mốc {i}: x = {x[i]}, y = {y[i]}")

    x0_str = input("\nMời nhập giá trị cần tính x0 (hoặc nhấn Enter để bỏ qua): ")
    x0 = float(x0_str) if x0_str.strip() != "" else None

    # Khởi tạo ma trận không 1 chiều
    BTH, BTT, f = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
    BTH[0], BTT[0], f[0] = y[0], 1, y[0]

    f, BTH, BTT = Newton_interpolation_Forward(BTH, BTT, f, x, y, n)

    value = fx0(f, n, x0) if x0 is not None else None

    print("\n--- KẾT QUẢ ĐA THỨC NỘI SUY NEWTON ---")
    dd_table = ComputeDividedDifferences(x, y)
    print(f"Dạng Newton: P(x) = {format_newton_form(x, dd_table[0, :])}")
    print(f"Dạng thu gọn: P(x) = {format_standard_poly(f)}")
    print(f"Hệ số từ x^0 đến x^{n}: {f}")

    if x0 is not None:
        print(f"Giá trị tính tại x0 = {x0} là: {value:.6f}")

    # Xuất ra file outputnewton.txt và vẽ hình
    ExportOutput(x, y, f, x0, value, is_added=False)
    plot(x, y, f, n)

    # Thêm mốc nội suy (nếu muốn)
    try:
        k_input = input("\nNhập số lượng thêm mốc nội suy (nhấn Enter để bỏ qua): ")
        if k_input.strip() == "":
            return
        k = int(k_input)
        if k <= 0:
            return
        for i in range(k):
            xt, yt = input(
                f"Nhập mốc và giá trị thứ {i+1} (cách nhau khoảng trắng): "
            ).split()
            if float(xt) not in x:
                n += 1
                x.append(float(xt))
                y.append(float(yt))
                f, BTH, BTT = add_Newton(BTH, BTT, f, x, y, n)

        value = fx0(f, n, x0) if x0 is not None else None
        print("\n--- ĐA THỨC MỚI SAU KHI THÊM MỐC ---")
        dd_table = ComputeDividedDifferences(x, y)
        print(f"Dạng Newton: P(x) = {format_newton_form(x, dd_table[0, :])}")
        print(f"Dạng thu gọn: P(x) = {format_standard_poly(f)}")

        if x0 is not None:
            print(f"Giá trị tính tại x0 = {x0} là: {value:.6f}")

        ExportOutput(x, y, f, x0, value, is_added=True)
        plot(x, y, f, n)
    except Exception as e:
        print("Kết thúc chương trình.")


if __name__ == "__main__":
    main()