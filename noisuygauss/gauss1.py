import math
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


# =====================================================================
# CÁC HÀM BỔ TRỢ ĐẠI SỐ ĐA THỨC & SAI PHÂN
# =====================================================================
def bangsaiphan(x, y):
    """
    Tính bảng sai phân: table[j][i] là sai phân cấp j tại mốc i.
    table[0] là các giá trị y ban đầu.
    """
    n = len(x)
    table = [[0.0] * n for _ in range(n)]
    for i in range(n):
        table[0][i] = float(y[i])

    for j in range(1, n):
        for i in range(n - j):
            table[j][i] = table[j - 1][i + 1] - table[j - 1][i]
    return table


def mulConst(poly, c):
    """Nhân đa thức với một hằng số c"""
    return [coef * c for coef in poly]


def sumHorner(poly1, poly2):
    """Cộng hai đa thức biểu diễn theo lũy thừa tăng dần của t (t^0, t^1, ...)"""
    max_len = max(len(poly1), len(poly2))
    res = [0.0] * max_len
    for i in range(len(poly1)):
        res[i] += poly1[i]
    for i in range(len(poly2)):
        res[i] += poly2[i]
    return res


def mulHorner(poly, step, direction=1):
    """
    Nhân đa thức tích của Gauss 1 theo từng bước:
    - Bước 1 (i = 1): nhân với t (tương đương t - 0)
    - Bước 2 (i = 2): nhân với (t - 1)
    - Bước 3 (i = 3): nhân với (t + 1)
    - Bước 4 (i = 4): nhân với (t - 2)
    - Bước 5 (i = 5): nhân với (t + 2)...
    """
    if step == 1:
        c = 0
    elif step % 2 == 0:
        c = step // 2
    else:
        c = -(step // 2)

    if len(poly) == 1 and poly[0] == 0:
        return [-c, 1.0]

    res = [0.0] * (len(poly) + 1)
    for i in range(len(poly)):
        res[i] -= c * poly[i]
        res[i + 1] += poly[i]
    return res


def horner(poly, t):
    """Tính giá trị đa thức P(t) bằng lược đồ Horner"""
    val = poly[-1]
    for i in range(len(poly) - 2, -1, -1):
        val = val * t + poly[i]
    return val


def format_poly_t(poly):
    """Định dạng đa thức thu gọn theo t: a_0 + a_1*t + ..."""
    terms = []
    for p, coef in enumerate(poly):
        if abs(coef) < 1e-9:
            continue
        val = abs(coef)
        sign = " - " if coef < 0 else " + "
        val_str = f"{val:g}" if abs(val - round(val)) < 1e-5 else f"{val:.5f}".rstrip("0").rstrip(".")

        if p == 0:
            term = f"{val_str}"
        elif p == 1:
            term = "t" if val == 1 else f"{val_str}*t"
        else:
            term = f"t^{p}" if val == 1 else f"{val_str}*t^{p}"

        if not terms:
            terms.append(f"-{term}" if coef < 0 else term)
        else:
            terms.append(f"{sign}{term}")
    return "".join(terms) if terms else "0"


def format_gauss1_factor_expr(list_sai_phan):
    """Định dạng biểu thức Gauss 1 dạng tích: y_0 + t*Δy_0 + t(t-1)/2!*Δ^2y_{-1} + ..."""
    terms = []
    fact = 1
    for i, sp in enumerate(list_sai_phan):
        if abs(sp) < 1e-9 and i > 0:
            continue
        if i > 1:
            fact *= i

        hs = sp / fact
        abs_hs = abs(hs)
        hs_str = f"{abs_hs:g}" if abs(abs_hs - round(abs_hs)) < 1e-5 else f"{abs_hs:.5f}".rstrip("0").rstrip(".")

        if i == 0:
            terms.append(f"{sp:.5f}".rstrip("0").rstrip("."))
        elif i == 1:
            sign = " + " if hs > 0 else " - "
            terms.append(f"{sign}{'t' if abs_hs == 1 else hs_str + '*t'}")
        else:
            factors = []
            for step in range(1, i + 1):
                if step == 1:
                    factors.append("t")
                elif step % 2 == 0:
                    factors.append(f"(t - {step // 2})")
                else:
                    factors.append(f"(t + {step // 2})")
            factor_str = "*".join(factors)
            sign = " + " if hs > 0 else " - "
            terms.append(f"{sign}{factor_str if abs_hs == 1 else hs_str + '*' + factor_str}")

    return "".join(terms) if terms else "0"


# =====================================================================
# THUẬT TOÁN GAUSS 1 VÀ XUẤT FILE BÁO CÁO BƯỚC LÀM
# =====================================================================
def gauss1_solver_and_export(x, y, value):
    n = len(x)
    h = x[1] - x[0]
    table = bangsaiphan(x, y)

    # Chọn mốc trung tâm x_0
    center_idx = int(n / 2) if (n % 2 == 1) else int((n - 1) / 2)
    x0_center = x[center_idx]
    t_val = (value - x0_center) / h

    list_sai_phan = []

    # Thuật toán tính hệ số đa thức Gauss 1
    poly = [0]
    ans = [table[0][center_idx]]
    fact = 1
    list_sai_phan.append(table[0][center_idx])

    for i in range(1, n):
        fact = fact * i
        poly = mulHorner(poly, i, 1)

        idx_row = int((n - i) / 2) if (n % 2 == 1) else int((n - i - 1) / 2)
        sp_val = table[i][idx_row]
        list_sai_phan.append(sp_val)

        term_poly = mulConst(poly, sp_val / fact)
        ans = sumHorner(ans, term_poly)

    result = horner(ans, t_val)

    # -------------------------------------------------------------
    # GHI CÁC BƯỚC LÀM VÀ KẾT QUẢ RA FILE outputgauss1.txt
    # -------------------------------------------------------------
    current_dir = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(current_dir, "outputgauss1.txt")

    lines = []
    lines.append("=" * 80)
    lines.append("               THUẬT TOÁN NỘI SUY GAUSS 1 (GAUSS THUẬN)")
    lines.append("=" * 80 + "\n")

    # 1. Bảng mốc
    lines.append(f"1. CÁC MỐC NỘI SUY (Bước nhảy h = {h:.6f}, Tổng số mốc = {n}):")
    lines.append(f"{'Chỉ số':<10}{'x_i':<15}{'y_i':<15}{'Ghi chú':<15}")
    lines.append("-" * 55)
    for i in range(n):
        note = "<- MỐC GỐC (x_0)" if i == center_idx else ""
        lines.append(f"x_{i - center_idx:<8}{x[i]:<15.5f}{y[i]:<15.5f}{note:<15}")
    lines.append("")

    # 2. Đổi biến t
    lines.append("2. CÔNG THỨC ĐỔI BIẾN SANG BIẾN t:")
    lines.append(f"   t = (x - x_0) / h = (x - {x0_center}) / {h:.6f}")
    lines.append(f"   Tại điểm cần tính x* = {value}:")
    lines.append(f"   => t* = ({value} - {x0_center}) / {h:.6f} = {t_val:.6f}\n")

    # 3. Bảng sai phân
    lines.append("3. BẢNG SAI PHÂN ĐẦY ĐỦ:")
    header = f"{'x_i':<10}{'y_i':<12}" + "".join([f"{f'Δ^{k}y':<12}" for k in range(1, n)])
    lines.append(header)
    lines.append("-" * len(header))
    for i in range(n):
        row = f"{x[i]:<10.4f}"
        for j in range(n - i):
            row += f"{table[j][i]:<12.4f}"
        lines.append(row)
    lines.append("-" * len(header) + "\n")

    # 4. Các sai phân được chọn theo đường zig-zag Gauss 1
    lines.append("4. CÁC SAI PHÂN ĐƯỢC CHỌN THEO QUY TẮC GAUSS 1 (ZIG-ZAG BẮT ĐẦU ĐI XUỐNG):")
    lines.append("   Đường đi: y_0 -> Δy_0 -> Δ^2 y_{-1} -> Δ^3 y_{-1} -> Δ^4 y_{-2} -> Δ^5 y_{-2} ...")
    for step, sp in enumerate(list_sai_phan):
        if step == 0:
            lines.append(f"   • Cấp 0: y_0 = {sp:.6f}")
        else:
            k = -(step // 2) if step % 2 == 0 else -( (step - 1) // 2 )
            lines.append(f"   • Cấp {step}: Δ^{step} y_{{{k}}} = {sp:.6f}")
    lines.append("")

    # 5. Biểu thức đa thức nội suy Gauss 1
    lines.append("5. ĐA THỨC NỘI SUY GAUSS 1:")
    lines.append("   Công thức tổng quát:")
    lines.append("   P(t) = y_0 + t*Δy_0 + [t(t-1)/2!]*Δ^2y_{-1} + [(t+1)t(t-1)/3!]*Δ^3y_{-1} + ...\n")
    lines.append("   a) Dạng tích thừa số Gauss 1 theo biến t:")
    lines.append(f"      P(t) = {format_gauss1_factor_expr(list_sai_phan)}\n")

    lines.append("   b) Dạng khai triển thu gọn theo lũy thừa của t:")
    lines.append(f"      P(t) = {format_poly_t(ans)}\n")

    # 6. Tính giá trị chi tiết tại value
    lines.append(f"6. TÍNH TOÁN GIÁ TRỊ TẠI x* = {value} (ứng với t* = {t_val:.6f}):")
    total_val = 0.0
    fact_k = 1
    term_mult = 1.0
    for k, sp in enumerate(list_sai_phan):
        if k == 0:
            term_val = sp
            lines.append(f"   • k = 0: y_0 = {term_val:.6f}")
        else:
            fact_k *= k
            if k == 1:
                term_mult = t_val
            elif k % 2 == 0:
                term_mult *= (t_val - (k // 2))
            else:
                term_mult *= (t_val + (k // 2))
            term_val = (sp / fact_k) * term_mult
            lines.append(f"   • k = {k}: (Δ^{k} y / {k}!) * [Tích thừa số t] = ({sp:.4f}/{fact_k}) * ({term_mult:.6f}) = {term_val:.6f}")
        total_val += term_val

    lines.append(f"\n=> KẾT QUẢ CUỐI CÙNG: P({value}) = {result:.10f}\n")

    # 7. Kiểm tra lại trên các mốc (Verification)
    lines.append("7. KIỂM TRA ĐỘ CHÍNH XÁC TẠI CÁC MỐC:")
    lines.append(f"{'x_i':<12}{'t_i':<12}{'y_thực tế':<15}{'y_nội suy':<15}{'Sai số':<15}")
    lines.append("-" * 69)
    for ele_x, ele_y in zip(x, y):
        ele_t = (ele_x - x0_center) / h
        pred = horner(ans, ele_t)
        err = abs(pred - ele_y)
        lines.append(f"{ele_x:<12.4f}{ele_t:<12.4f}{ele_y:<15.5f}{pred:<15.5f}{err:<15.2e}")

    content = "\n".join(lines)
    with open(out_file, "w", encoding="utf-8") as f_out:
        f_out.write(content + "\n")

    print(content)
    print(f"\n>>> Đã xuất toàn bộ kết quả và các bước giải vào: {out_file}")
    return ans, result


# =====================================================================
# HÀM CHÍNH (MAIN)
# =====================================================================
if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "vd1.txt")

    x = []
    y = []

    # 1. Đọc dữ liệu từ file vd1.txt
    if os.path.exists(input_file):
        with open(input_file, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
    else:
        print(f"Không tìm thấy file {input_file}, dùng dữ liệu mẫu mặc định:")
        x = [310, 320, 330, 340, 350, 360, 370]
        y = [2.4914, 2.5052, 2.5185, 2.5315, 2.5441, 2.5563, 2.5687]

    print("Dữ liệu các mốc đã đọc:")
    for xi, yi in zip(x, y):
        print(f"  x = {xi:<8.2f} y = {yi:<8.4f}")

    # 2. Nhập giá trị cần tính từ bàn phím (có giá trị mặc định 345 nếu nhấn Enter)
    user_input = input("\nNhập giá trị cần tính nội suy (nhấn Enter để dùng mặc định 345): ").strip()
    if user_input == "":
        value = 345.0
    else:
        value = float(user_input)

    print(f"\n>>> Đang tính nội suy tại x = {value}...")

    # 3. Chạy thuật toán và xuất báo cáo ra outputgauss1.txt
    ans, result = gauss1_solver_and_export(x, y, value)

    # 4. Vẽ đồ thị minh họa
    n = len(x)
    center_idx = int(n / 2) if (n % 2 == 1) else int((n - 1) / 2)
    x0_center = x[center_idx]
    h = x[1] - x[0]

    plt.figure(figsize=(9, 5))
    plt.plot(x, y, "ro", label="Dữ liệu mốc (Data points)")

    x_plot = np.linspace(min(x), max(x), 200)
    t_plot = (x_plot - x0_center) / h
    y_plot = [horner(ans, t) for t in t_plot]

    plt.plot(x_plot, y_plot, "b-", label="Đa thức nội suy Gauss 1")
    plt.plot([value], [result], "g*", markersize=12, label=f"P({value}) = {result:.5f}")

    plt.title("Nội suy Gauss 1 (Gauss Forward)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()

    img_file = os.path.join(current_dir, "gauss1.png")
    plt.savefig(img_file, dpi=300)
    plt.show()
    print(f">>> Đã lưu đồ thị vào: {img_file}")