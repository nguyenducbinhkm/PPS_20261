import os
import sys


class Data:

    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y


def format_poly(coeffs):
    """
    Hàm định dạng in đa thức toán học đẹp mắt.
    coeffs: danh sách hệ số từ x^0 đến x^m (chỉ số i là bậc x^i)
    """
    terms = []
    for p, coef in enumerate(coeffs):
        if abs(coef) < 1e-9:
            continue
        val = abs(coef)
        sign = " - " if coef < 0 else " + "
        val_str = f"{val:g}" if abs(val - round(val)) < 1e-6 else f"{val:.4f}"

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


def hoocnerNhan(arr, hs, xk):
    """Nhân đa thức arr với nhị thức (hs * x + xk)"""
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


def hoocnerChia(arr, hs, xk):
    """Chia đa thức arr cho nhị thức (hs * x + xk) bằng sơ đồ Horner"""
    ret_size = len(arr) - 1
    ret = [0.0] * ret_size
    nho = 0.0
    for i in range(len(arr) - 1, 0, -1):
        val = (arr[i] - nho) / hs
        nho = val * xk
        ret[i - 1] = val
    return ret


class DualLogger:
    """Ghi đồng thời ra Terminal và File out"""

    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.log_file = open(filepath, "w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)

    def flush(self):
        self.terminal.flush()
        self.log_file.flush()


def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    input_path = os.path.join(current_dir, "nsLagrange.inp")
    output_path = os.path.join(current_dir, "nsLagrange.out")

    try:
        with open(input_path, "r", encoding="utf-8") as f:
            input_data = f.read().split()
    except FileNotFoundError:
        print(f"Lỗi: Không tìm thấy file {input_path}")
        return

    if not input_data:
        return

    # Kích hoạt bộ ghi kép (console + file out.txt)
    logger = DualLogger(output_path)
    sys.stdout = logger

    it = iter(input_data)
    num_points = int(next(it))  # Số lượng mốc nội suy
    n = num_points - 1  # Bậc tối đa của đa thức nội suy

    a = [Data() for _ in range(num_points)]
    for i in range(num_points):
        a[i].x = float(next(it))
        a[i].y = float(next(it))

    print("=" * 75)
    print("      THUẬT TOÁN XÂY DỰNG ĐA THỨC NỘI SUY LAGRANGE CHI TIẾT")
    print("=" * 75)

    # -------------------------------------------------------------
    # BƯỚC 1: SẮP XẾP CÁC MỐC NỘI SUY
    # -------------------------------------------------------------
    a.sort(key=lambda item: item.x)
    print("\n[BƯỚC 1] BẢNG CÁC MỐC NỘI SUY (ĐÃ SẮP XẾP TĂNG DẦN THEO X):")
    print(f"{'Mốc i':<8}{'x_i':<15}{'y_i = f(x_i)':<15}")
    print("-" * 38)
    for i in range(num_points):
        print(f"{i:<8}{a[i].x:<15.4f}{a[i].y:<15.4f}")
    print("-" * 38)

    # -------------------------------------------------------------
    # BƯỚC 2: XÂY DỰNG ĐA THỨC TÍCH w(x) = (x - x_0)...(x - x_n)
    # -------------------------------------------------------------
    print("\n[BƯỚC 2] XÂY DỰNG ĐA THỨC TÍCH:")
    print("   ω(x) = (x - x_0)(x - x_1)...(x - x_n)")

    arr_wi = []
    for i in range(num_points):
        arr_wi = hoocnerNhan(arr_wi, 1.0, -a[i].x)
        print(f"   • Nhân thêm (x - {a[i].x}):")
        print(f"     => ω_{i}(x) = {format_poly(arr_wi)}")

    print(f"\n=> ĐA THỨC TÍCH TOÀN PHẦN:")
    print(f"   ω(x) = {format_poly(arr_wi)}\n")

    # -------------------------------------------------------------
    # BƯỚC 3: TÍNH ĐẠO HÀM w'(x_i) = TÍCH (x_i - x_j) với j != i
    # -------------------------------------------------------------
    print("\n[BƯỚC 3] TÍNH ĐẠO HÀM ω'(x_i) TẠI TỪNG MỐC NỘI SUY:")
    print("   Công thức: ω'(x_i) = ∏_{j ≠ i} (x_i - x_j)\n")

    dw = [1.0] * num_points
    for i in range(num_points):
        step_str = []
        for j in range(num_points):
            if i != j:
                diff = a[i].x - a[j].x
                dw[i] *= diff
                step_str.append(f"({a[i].x} - {a[j].x})")
        print(f"   • Tại x_{i} = {a[i].x}:")
        print(f"     ω'({a[i].x}) = {' * '.join(step_str)} = {dw[i]:.6f}")

    # -------------------------------------------------------------
    # BƯỚC 4: TÍNH CÁC ĐA THỨC CƠ SỞ L_i(x) VÀ HỆ SỐ ĐÓNG GÓP
    # -------------------------------------------------------------
    print("\n[BƯỚC 4] XÁC ĐỊNH CÁC ĐA THỨC CƠ SỞ VÀ THÀNH PHẦN:")
    print("   Công thức mốc cơ sở: l_i(x) = ω(x) / [(x - x_i) * ω'(x_i)]")
    print("   Thành phần đóng góp: P_i(x) = y_i * l_i(x) = [y_i / ω'(x_i)] * [ω(x) / (x - x_i)]\n")

    arr_Lagrange = [0.0] * num_points
    for i in range(num_points):
        # Chia w(x) cho (x - x_i) bằng sơ đồ Horner
        wi_div = hoocnerChia(arr_wi, 1.0, -a[i].x)
        hs = a[i].y / dw[i]

        print(f"   --- [Mốc i = {i}: x_{i} = {a[i].x}, y_{i} = {a[i].y}] ---")
        print(f"   • ω(x)/(x - {a[i].x}) = {format_poly(wi_div)}")
        print(f"   • Tỷ số y_{i} / ω'({a[i].x}) = {a[i].y} / {dw[i]:.6f} = {hs:.6f}")
        print(f"   • Đóng góp P_{i}(x) = {hs:.6f} * [{format_poly(wi_div)}]")

        # Cộng dồn vào đa thức Lagrange tổng
        for k in range(num_points):
            arr_Lagrange[k] += wi_div[k] * hs
        print()

    # -------------------------------------------------------------
    # BƯỚC 5: TỔNG HỢP ĐA THỨC NỘI SUY LAGRANGE
    # -------------------------------------------------------------
    print("=" * 75)
    print("[BƯỚC 5] ĐA THỨC NỘI SUY LAGRANGE TỔNG HỢP: L_n(x) = ∑ P_i(x)")
    print("=" * 75)
    print("Các hệ số thu gọn từ bậc x^0 đến bậc x^n:")
    for deg, val in enumerate(arr_Lagrange):
        print(f"   • Hệ số của x^{deg}: {val:.8f}")

    print(f"\n=> DẠNG BIỂU THỨC ĐẠI SỐ ĐA THỨC NỘI SUY:")
    print(f"   P(x) = {format_poly(arr_Lagrange)}\n")

    # -------------------------------------------------------------
    # BƯỚC 6: TÍNH GIÁ TRỊ TẠI CÁC ĐIỂM TRUY VẤN (NẾU CÓ)
    # -------------------------------------------------------------
    try:
        nq_token = next(it, None)
        if nq_token is not None:
            nq = int(nq_token)
            print("=" * 75)
            print(f"[BƯỚC 6] TÍNH GIÁ TRỊ TẠI {nq} ĐIỂM TRUY VẤN BẰNG P(x):")
            print("=" * 75)

            for q_idx in range(1, nq + 1):
                value_x = float(next(it))
                # Thay số tính từng bước: sum(a_k * x^k)
                val_y = 0.0
                step_calc = []
                for k in range(len(arr_Lagrange)):
                    term_val = arr_Lagrange[k] * (value_x**k)
                    val_y += term_val
                    step_calc.append(f"({arr_Lagrange[k]:.4f} * {value_x}^{k})")

                print(f"\n• Truy vấn {q_idx}: Tính giá trị tại x = {value_x}")
                print(f"  P({value_x}) = " + " + ".join(step_calc))
                print(f"  => P({value_x}) = {val_y:.10f}")
    except StopIteration:
        pass

    print("\n" + "=" * 75)
    print(f">>> HOÀN TẤT! Toàn bộ các bước trên đã được ghi vào: nsLagrange.out")
    print("=" * 75)


if __name__ == "__main__":
    main()