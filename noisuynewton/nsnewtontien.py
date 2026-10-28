from math import *
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(suppress=True, linewidth=np.inf, precision=5)


def Input():
    x, y = [], []
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, "input.txt")

    if not os.path.exists(file_path):
        print(f"Lỗi: Không tìm thấy file {file_path}")
        sys.exit()

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f.readlines():
            if not line.strip():
                continue
            tokens = line.split()
            val_x = float(tokens[0])
            if val_x not in x:
                x.append(val_x)
                y.append(float(tokens[1]))

    # Kiểm tra điều kiện mốc cách đều
    if len(x) < 2:
        print("Lỗi: Cần ít nhất 2 mốc nội suy trong file input.txt.")
        sys.exit()

    h = x[1] - x[0]
    for i in range(1, len(x)):
        if abs(x[i] - x[i - 1] - h) > 1e-6:
            print("Input không hợp lệ! Các mốc phải cách đều nhau.")
            sys.exit()

    return x, y


# Chọn ra 'num' mốc gần x0 nhất
def pickPoint(x0, num, x, y):
    if num > len(x):
        print(f"Lỗi: Số lượng mốc yêu cầu ({num}) lớn hơn tổng số mốc trong file ({len(x)}).")
        sys.exit()
    h = x[1] - x[0]
    i = int((x0 - x[0]) / h)
    left = min(len(x) - num, max(0, i + 1 - int(num / 2)))
    right = left + num - 1
    return x[left : right + 1], y[left : right + 1]


# Tính bảng sai phân tiến đầy đủ
def ComputeForwardDifferenceTable(y):
    n = len(y)
    diff_table = np.zeros((n, n))
    for i in range(n):
        diff_table[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            diff_table[i][j] = diff_table[i + 1][j - 1] - diff_table[i][j - 1]

    return diff_table


# Kết nạp mốc x[i] vào Bảng sai phân Tiến (theo code gốc)
def BuildBSP(BSP, x, y, i):
    SP_row_i = np.zeros(len(x))
    SP_row_i[0] = y[i]
    for j in range(1, i + 1):
        SP_row_i[j] = SP_row_i[j - 1] - BSP[j - 1]
    return SP_row_i


# Kết nạp mốc x[i] vào Bảng tính tích biến t
def BuildBTT(BTT, x, i):
    TT_row_i = np.zeros(len(x))
    TT_row_i[i] = 1
    for j in range(i, 0, -1):
        TT_row_i[j] = BTT[j - 1] - i * BTT[j]
    if i != 0:
        TT_row_i[0] = -i * BTT[0]
    return TT_row_i


# Kết nạp mốc x[i] vào đa thức Newton Tiến
def fx(BSP, BTT, i, f):
    for j in range(1, i + 1):
        f[j] += BSP[i] * BTT[j - 1] / factorial(i)
    return f


# Xây dựng đa thức nội suy Newton tiến
def Newton_interpolation_Forward(BSP, BTT, f, x1, y1, n):
    for i in range(1, n + 1):
        BSP = BuildBSP(BSP, x1, y1, i)
        BTT = BuildBTT(BTT, x1, i - 1)
        f = fx(BSP, BTT, i, f)
    return f, BSP, BTT


# Ước tính f(x0)
def fx0(f, x0, x1):
    t = (x0 - x1[0]) / (x1[1] - x1[0])
    value = f[len(x1) - 1]
    for i in range(len(x1) - 2, -1, -1):
        value = value * t + f[i]
    return value


# Định dạng đa thức theo biến t
def format_poly_t(coeffs):
    terms = []
    for p, c in enumerate(coeffs):
        if abs(c) < 1e-9:
            continue
        val = abs(c)
        sign = " - " if c < 0 else " + "
        val_str = f"{val:g}" if abs(val - round(val)) < 1e-5 else f"{val:.5f}".rstrip("0").rstrip(".")

        if p == 0:
            term = f"{val_str}"
        elif p == 1:
            term = "t" if val == 1 else f"{val_str}*t"
        else:
            term = f"t^{p}" if val == 1 else f"{val_str}*t^{p}"

        if not terms:
            terms.append(f"-{term}" if c < 0 else term)
        else:
            terms.append(f"{sign}{term}")

    return "".join(terms) if terms else "0"


# Định dạng tích Newton tiến: P(t) = y0 + dy0*t + d2y0/2!*t(t-1) + ...
def format_newton_forward_expr(delta_y):
    terms = []
    for k, val in enumerate(delta_y):
        if abs(val) < 1e-9 and k > 0:
            continue
        hs = val / factorial(k)
        if k == 0:
            terms.append(f"{val:.5f}".rstrip("0").rstrip("."))
        else:
            factors = "*".join([f"(t - {j})" if j > 0 else "t" for j in range(k)])
            sign = " + " if hs > 0 else " - "
            abs_hs = abs(hs)
            hs_str = f"{abs_hs:g}" if abs(abs_hs - round(abs_hs)) < 1e-5 else f"{abs_hs:.5f}".rstrip("0").rstrip(".")
            term = factors if abs_hs == 1 else f"{hs_str}*{factors}"
            terms.append(f"{sign}{term}")

    return "".join(terms) if terms else "0"


# Hàm xuất báo cáo chi tiết thuật toán ra outputnewtontien.txt
def ExportOutput(x1, y1, f, x0, is_added=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(current_dir, "outputnewtontien.txt")
    mode = "a" if is_added else "w"

    h = x1[1] - x1[0]
    t0 = (x0 - x1[0]) / h
    n = len(x1) - 1
    diff_table = ComputeForwardDifferenceTable(y1)
    delta_y0 = [diff_table[0][k] for k in range(n + 1)]

    lines = []
    if not is_added:
        lines.append("=" * 80)
        lines.append("           THUẬT TOÁN NỘI SUY NEWTON TIẾN (MỐC CÁCH ĐỀU)")
        lines.append("=" * 80 + "\n")
    else:
        lines.append("\n" + "=" * 80)
        lines.append("       KẾT QUẢ NỘI SUY NEWTON TIẾN SAU KHI THÊM MỐC MỚI")
        lines.append("=" * 80 + "\n")

    # 1. Bảng mốc
    lines.append(f"1. CÁC MỐC NỘI SUY ĐÃ CHỌN (Bước nhảy h = {h:.6f}):")
    lines.append(f"{'i':<5}{'x_i':<15}{'y_i':<15}")
    lines.append("-" * 35)
    for i in range(len(x1)):
        lines.append(f"{i:<5}{x1[i]:<15.5f}{y1[i]:<15.5f}")
    lines.append("")

    # 2. Đổi biến t
    lines.append("2. CÔNG THỨC ĐỔI BIẾN KHÔNG THỨ NGUYÊN:")
    lines.append(f"   t = (x - x_0) / h = (x - {x1[0]}) / {h:.6f}")
    lines.append(f"   Tại điểm cần tính x* = {x0}:")
    lines.append(f"   => t* = ({x0} - {x1[0]}) / {h:.6f} = {t0:.6f}\n")

    # 3. Bảng sai phân tiến
    lines.append("3. BẢNG SAI PHÂN TIẾN ĐẦY ĐỦ:")
    header = f"{'i':<4}{'x_i':<10}{'y_i':<12}" + "".join([f"{f'Δ^{k}y':<12}" for k in range(1, n + 1)])
    lines.append(header)
    lines.append("-" * len(header))
    for i in range(len(x1)):
        row = f"{i:<4}{x1[i]:<10.4f}{y1[i]:<12.4f}"
        for j in range(1, len(x1) - i):
            row += f"{diff_table[i][j]:<12.4f}"
        lines.append(row)
    lines.append("-" * len(header))

    lines.append("\n   Các sai phân tại mốc ban đầu (hàng 0) dùng cho công thức Newton tiến:")
    for k in range(n + 1):
        lines.append(f"   • Δ^{k} y_0 = {delta_y0[k]:.6f}")
    lines.append("")

    # 4. Biểu thức đa thức
    lines.append("4. BIỂU THỨC ĐA THỨC NỘI SUY NEWTON TIẾN:")
    lines.append("   Công thức tổng quát:")
    lines.append("   P(t) = y_0 + (Δy_0/1!)*t + (Δ^2 y_0/2!)*t(t-1) + ... + (Δ^n y_0/n!)*t(t-1)...(t-n+1)\n")
    lines.append(f"   a) Dạng tích Newton tiến theo biến t:")
    lines.append(f"      P(t) = {format_newton_forward_expr(delta_y0)}\n")

    lines.append(f"   b) Dạng khai triển thu gọn theo lũy thừa của t:")
    lines.append(f"      P(t) = {format_poly_t(f)}\n")

    # 5. Các bước thay số tính f(x0)
    lines.append(f"5. TÍNH TOÁN GIÁ TRỊ TẠI ĐIỂM x* = {x0} (ứng với t* = {t0:.6f}):")
    total_val = 0.0
    t_prod = 1.0
    lines.append("   Chi tiết từng số hạng trong công thức Newton:")
    for k in range(n + 1):
        if k == 0:
            term_val = delta_y0[0]
            lines.append(f"   • k = 0: y_0 = {term_val:.6f}")
        else:
            t_prod *= (t0 - (k - 1))
            term_val = (delta_y0[k] / factorial(k)) * t_prod
            lines.append(f"   • k = {k}: (Δ^{k}y_0 / {k}!) * [t*(t-1)...(t-{k-1})] = ({delta_y0[k]:.4f}/{factorial(k)}) * ({t_prod:.6f}) = {term_val:.6f}")
        total_val += term_val

    lines.append(f"\n=> KẾT QUẢ CUỐI CÙNG: P({x0}) = {total_val:.10f}")

    content = "\n".join(lines)
    with open(out_file, mode, encoding="utf-8") as f_out:
        f_out.write(content + "\n")

    print(content)
    print(f"\n>>> Đã lưu toàn bộ các bước làm vào file: {out_file}")


# Thêm mốc nội suy
def add_Newton(BSP, BTT, f, x1, y1, n):
    BSP = BuildBSP(BSP, x1, y1, n)
    BTT = BuildBTT(BTT, x1, n - 1)
    f = np.hstack((f, np.zeros(1)))
    f = fx(BSP, BTT, n, f)
    return f, BSP, BTT


# Vẽ hình
def plot(x, y, f):
    plt.figure()
    plt.title("Newton Forward Interpolation")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.scatter(x, y, color="red", label="Data points")
    xx = np.linspace(min(x), max(x), 1000)
    plt.plot(xx, [fx0(f, val, x) for val in xx], color="blue", label="P_Newton_Forward(x)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    current_dir = os.path.dirname(os.path.abspath(__file__))
    img_path = os.path.join(current_dir, "newton_forward.png")
    plt.savefig(img_path, dpi=300)
    plt.close()
    print(f">>> Đã lưu đồ thị vào file: {img_path}")


def main():
    x, y = Input()
    x0 = float(input("Mời nhập giá trị cần tính x0: "))
    num = int(input(f"Chọn số lượng mốc nội suy tính (<= {len(x)}): "))

    h = x[1] - x[0]
    x1, y1 = pickPoint(x0, num, x, y)
    n = len(x1) - 1

    BSP, BTT, f = np.zeros(n + 1), np.zeros(n + 1), np.zeros(n + 1)
    BSP[0], BTT[0], f[0] = y1[0], 1, y1[0]

    f, BSP, BTT = Newton_interpolation_Forward(BSP, BTT, f, x1, y1, n)

    # Xuất các bước làm và kết quả ra file outputnewtontien.txt
    ExportOutput(x1, y1, f, x0, is_added=False)
    plot(x1, y1, f)

    try:
        k_input = input("\nNhập số lượng thêm mốc nội suy (nhấn Enter hoặc 0 để bỏ qua): ")
        if k_input.strip() == "" or int(k_input) <= 0:
            return
        k = int(k_input)
        for i in range(k):
            xt, yt = input(f"Nhập mốc và giá trị thứ {i+1}: ").split()
            if float(xt) not in x1 and abs(float(xt) - x1[n] - h) < 1e-6:
                n += 1
                x1.append(float(xt))
                y1.append(float(yt))
                f, BSP, BTT = add_Newton(BSP, BTT, f, x1, y1, n)

        ExportOutput(x1, y1, f, x0, is_added=True)
        plot(x1, y1, f)
    except Exception:
        print("Kết thúc phần thêm mốc.")


if __name__ == "__main__":
    main()