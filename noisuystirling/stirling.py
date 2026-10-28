from math import factorial, prod
import os
import sys
import numpy as np

np.set_printoptions(precision=10, suppress=True)


# =====================================================================
# CÁC HÀM BỔ TRỢ ĐẠI SỐ ĐA THỨC THEO CÔNG THỨC VIETE
# =====================================================================
def vietaFormula(coeff: list, x0: float) -> list:
    """Nhân đa thức coeff với (t^2 - x0)"""
    coeff = [0] + coeff
    n = len(coeff)
    for i in range(n - 1):
        coeff[i] -= x0 * coeff[i + 1]
    return coeff


def plus_poly(high: list, low: list) -> list:
    """Cộng 2 đa thức hệ số"""
    n_high = len(high)
    n_low = len(low)
    res = [0.0] * max(n_high, n_low)
    for i in range(n_high):
        res[i] += high[i]
    for i in range(n_low):
        res[i] += low[i]
    return res


def mul_list(kirito: list, ratio: float) -> list:
    return [i * ratio for i in kirito]


def deri(ys: list):
    """Tính sai phân cấp 1 từ dãy ys"""
    n = len(ys)
    return [ys[i + 1] - ys[i] for i in range(n - 1)]


def format_poly_t(coeffs):
    """Định dạng đa thức theo lũy thừa tăng dần của t: a0 + a1*t + a2*t^2 + ..."""
    terms = []
    for p, c in enumerate(coeffs):
        if abs(c) < 1e-11:
            continue
        val = abs(c)
        sign = " - " if c < 0 else " + "
        val_str = (
            f"{val:g}"
            if abs(val - round(val)) < 1e-6
            else f"{val:.6f}".rstrip("0").rstrip(".")
        )

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


# =====================================================================
# CHỌN MỐC NỘI SUY TỐI ƯU ĐỐI XỨNG QUANH ĐIỂM CẦN TÍNH
# =====================================================================
def nodes(xs: list, ys: list, c: float):
    p = len(xs)
    h = xs[1] - xs[0]
    j = 0

    if abs(c - xs[0]) < abs(c - xs[-1]):
        i = 0
        while i < p and abs(xs[i] - c) > h / 2:
            i += 1
        if i < p and xs[i] == c:
            print(f"Giá trị tại {c} đã biết chính xác là: {ys[i]}")
            sys.exit()

        while (
            0 <= i + j <= p - 1
            and 0 <= i - j <= p - 1
            and xs[i + j] - xs[i - j] <= 1
        ):
            j += 1
        j -= 1
        if j <= 1:
            if i >= 2 and i + 2 < p:
                return range(i - 2, i + 3)
            else:
                print("Lỗi: Quá ít điểm dữ liệu xung quanh mốc cần tính để nội suy Stirling.")
                sys.exit()
        elif 2 <= j <= 4:
            return range(i - j, i + j + 1)
        else:
            return range(i - 4, i + 5)
    else:
        i = p - 1
        while i >= 0 and abs(xs[i] - c) > h / 2:
            i -= 1
        if i >= 0 and xs[i] == c:
            print(f"Giá trị tại {c} đã biết chính xác là: {ys[i]}")
            sys.exit()

        while (
            0 <= i + j <= p - 1
            and 0 <= i - j <= p - 1
            and xs[i + j] - xs[i - j] <= 1
        ):
            j += 1
        j -= 1
        if j <= 1:
            if i >= 2 and i + 2 < p:
                return range(i - 2, i + 3)
            else:
                print("Lỗi: Quá ít điểm dữ liệu xung quanh mốc cần tính để nội suy Stirling.")
                sys.exit()
        elif 2 <= j <= 4:
            return range(i - j, i + j + 1)
        else:
            return range(i - 4, i + 5)


# =====================================================================
# THUẬT TOÁN STIRLING OVERHAUL & XUẤT BÁO CÁO CHI TIẾT
# =====================================================================
def Stirling_overhaul(all_xs: list, all_ys: list, x_val: float, output):
    selected_indices = list(nodes(all_xs, all_ys, x_val))
    xs = [all_xs[idx] for idx in selected_indices]
    ys = [all_ys[idx] for idx in selected_indices]

    total_pts = len(xs)
    n = (total_pts - 1) // 2  # Bán kính mốc đối xứng quanh x0: x_{-n} ... x_0 ... x_n
    h = round(xs[1] - xs[0], 10)
    x0 = xs[n]
    y0 = ys[n]
    t = (x_val - x0) / h

    # Lập bảng sai phân đầy đủ để hiển thị
    diff_table = [ys[:]]
    curr_diff = ys[:]
    for _ in range(total_pts - 1):
        curr_diff = deri(curr_diff)
        diff_table.append(curr_diff)

    def pr(text=""):
        print(text, file=output)
        print(text)

    pr("=" * 80)
    pr(f"                 THUẬT TOÁN NỘI SUY STIRLING TẠI x = {x_val}")
    pr("=" * 80 + "\n")

    # 1. Các mốc được chọn
    pr(f"1. CÁC MỐC NỘI SUY ĐỐI XỨNG ĐÃ CHỌN (Tổng cộng 2n + 1 = {total_pts} mốc, Bước h = {h}):")
    pr(f"{'Chỉ số':<10}{'x_i':<15}{'y_i = f(x_i)':<18}{'Ghi chú':<18}")
    pr("-" * 60)
    for idx, (xi, yi) in enumerate(zip(xs, ys)):
        relative_idx = idx - n
        note = "<- MỐC GỐC TRUNG TÂM (x_0)" if relative_idx == 0 else ""
        pr(f"x_{relative_idx:<8}{xi:<15.6f}{yi:<18.8f}{note:<18}")
    pr("")

    # 2. Đổi biến không thứ nguyên t
    pr("2. CÔNG THỨC ĐỔI BIẾN SANG BIẾN t:")
    pr(f"   t = (x - x_0) / h = (x - {x0}) / {h}")
    pr(f"   Tại x* = {x_val}:")
    pr(f"   => t* = ({x_val} - {x0}) / {h} = {t:.10f}\n")

    # 3. Bảng sai phân trung tâm
    pr("3. BẢNG SAI PHÂN ĐẦY ĐỦ:")
    header = f"{'Mốc':<8}{'x_i':<10}{'y_i':<14}" + "".join(
        [f"{f'Δ^{k}y':<14}" for k in range(1, total_pts)]
    )
    pr(header)
    pr("-" * len(header))
    for i in range(total_pts):
        rel_idx = i - n
        row_str = f"x_{rel_idx:<6}{xs[i]:<10.4f}{ys[i]:<14.6f}"
        for col in range(1, total_pts):
            if i < len(diff_table[col]):
                row_str += f"{diff_table[col][i]:<14.6f}"
        pr(row_str)
    pr("-" * len(header) + "\n")

    # 4. Xác định các hệ số a_i và b_i của Stirling
    pr("4. CÁC HỆ SỐ SAI PHÂN TRUNG BÌNH STIRLING (a_k và b_k):")
    pr("   • Bậc lẻ  : a_k = [Δ^(2k-1) y_{-(k-1)} + Δ^(2k-1) y_{-k}] / [2 * (2k-1)!]")
    pr("   • Bậc chẵn: b_k = [Δ^(2k) y_{-k}] / [(2k)!]\n")

    coeff = [1]
    even_coeff = []
    odd_coeff = []
    fact = 1
    terms_calc_detail = []

    # Số hạng đầu tiên: y0
    terms_calc_detail.append((0, "y_0", y0))

    temp_ys = ys[:]
    for i in range(n):
        k = i + 1  # bước k = 1, 2, ..., n
        # Sai phân cấp lẻ: 2k - 1
        temp_ys = deri(temp_ys)
        fact_odd = fact * (2 * i + 1)
        # Các sai phân tương ứng
        sp_odd_1 = temp_ys[n - i]
        sp_odd_2 = temp_ys[n - i - 1]
        a_k = (sp_odd_1 + sp_odd_2) / (2 * fact_odd)

        pr(f"   --- BẬC {2*k - 1} (LẺ): k = {k} ---")
        pr(f"   • Cặp sai phân cấp {2*k - 1}: Δ^{2*k - 1}y_{{-(k-1)}} = {sp_odd_1:.8f}, Δ^{2*k - 1}y_{{-k}} = {sp_odd_2:.8f}")
        pr(f"   • a_{k} = ({sp_odd_1:.8f} + {sp_odd_2:.8f}) / [2 * {fact_odd}] = {a_k:.10f}")

        odd_coeff = plus_poly(mul_list(coeff, a_k), odd_coeff)

        # Sai phân cấp chẵn: 2k
        temp_ys = deri(temp_ys)
        fact_even = fact_odd * (2 * i + 2)
        fact = fact_even
        sp_even = temp_ys[n - i - 1]
        b_k = sp_even / fact_even

        pr(f"   --- BẬC {2*k} (CHẴN): k = {k} ---")
        pr(f"   • Sai phân cấp {2*k}: Δ^{2*k}y_{{-k}} = {sp_even:.8f}")
        pr(f"   • b_{k} = {sp_even:.8f} / {fact_even} = {b_k:.10f}\n")

        even_coeff = plus_poly(mul_list(coeff, b_k), even_coeff)
        coeff = vietaFormula(coeff, k**2)

    # 5. Xây dựng đa thức thu gọn theo t
    # Tổng hợp các hệ số đa thức P(t) = c0 + c1*t + c2*t^2 + ...
    final_poly = [y0]
    deg = 2 * n
    for j in range(n):
        c_odd = odd_coeff[j] if j < len(odd_coeff) else 0.0
        c_even = even_coeff[j] if j < len(even_coeff) else 0.0
        final_poly.append(c_odd)
        final_poly.append(c_even)

    pr("5. DẠNG BIỂU THỨC ĐA THỨC NỘI SUY STIRLING:")
    pr("   a) Công thức dạng tích Stirling:")
    pr("      P(t) = y_0 + a_1*t + b_1*t^2 + a_2*t(t^2 - 1^2) + b_2*t^2(t^2 - 1^2) + ...\n")

    pr("   b) Dạng khai triển thu gọn theo lũy thừa của t:")
    pr(f"      P(t) = {format_poly_t(final_poly)}\n")

    # 6. Tính giá trị chi tiết tại x*
    pr(f"6. CHI TIẾT TÍNH TOÁN GIÁ TRỊ TẠI x* = {x_val} (ứng với t* = {t:.8f}):")
    running_sum = 0.0
    for p, c in enumerate(final_poly):
        term_val = c * (t**p)
        running_sum += term_val
        pr(f"   • Số hạng bậc t^{p:<2}: {c:>15.10f} * ({t:.6f})^{p} = {term_val:>18.12f}")

    pr(f"\n=> KẾT QUẢ NỘI SUY STIRLING: P({x_val}) = {running_sum:.12f}\n")

    # 7. Ước lượng sai số
    last_diff = temp_ys[0]
    t_sq = t**2
    error_bound = abs(
        (last_diff / fact) * t * prod((t_sq - (k + 1) ** 2 for k in range(n)))
    )

    pr("7. ĐÁNH GIÁ SAI SỐ ƯỚC LƯỢNG:")
    pr("   Công thức: R_{2n}(x) ≈ [Δ^(2n+1) y / (2n+1)!] * t * ∏_{k=1}^n (t^2 - k^2)")
    pr(f"   => Sai số tuyệt đối |R(x)| không vượt quá: {error_bound:.10e}")
    pr("=" * 80)


# =====================================================================
# HÀM CHÍNH
# =====================================================================
def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_dir, "interpolation.txt")
    output_file = os.path.join(current_dir, "output.txt")

    if not os.path.exists(input_file):
        print(f"Lỗi: Không tìm thấy file {input_file}")
        return

    xs = []
    ys = []
    with open(input_file, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Tự động nhận diện dòng mốc dữ liệu và giá trị start cần tính
    # Dòng cuối là điểm start cần tính nếu chỉ có 1 số, hoặc các dòng có 2 số là mốc
    start_val = None
    last_tokens = lines[-1].split()
    if len(last_tokens) == 1:
        start_val = float(last_tokens[0])
        data_lines = lines[:-1]
    else:
        data_lines = lines

    for line in data_lines:
        tokens = line.replace(",", " ").split()
        if len(tokens) >= 2:
            xs.append(float(tokens[0]))
            ys.append(float(tokens[1]))

    if start_val is None:
        user_in = input(f"Nhập điểm cần tính nội suy x* (nằm trong khoảng [{xs[0]}, {xs[-1]}]): ").strip()
        start_val = float(user_in)

    if not (xs[0] <= start_val <= xs[-1]):
        print(f"Lỗi: Điểm cần tính {start_val} nằm ngoài khoảng [{xs[0]}, {xs[-1]}]. Phương pháp không hỗ trợ ngoại suy.")
        return

    with open(output_file, "w", encoding="utf-8") as output:
        Stirling_overhaul(xs, ys, start_val, output)

    print(f"\n>>> Đã xuất báo cáo các bước giải chi tiết vào file: {output_file}")


if __name__ == "__main__":
    main()