import math
import os
import sys

# Thiết lập đường dẫn file inp.txt và outdaoham.txt cùng thư mục với script
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
INP_PATH = os.path.join(CURRENT_DIR, "inp.txt")
OUT_PATH = os.path.join(CURRENT_DIR, "outdaoham.txt")


def format_poly(coeffs):
    """Hàm hỗ trợ in đa thức đẹp mắt từ danh sách hệ số [a_m, a_{m-1}, ..., a_0]"""
    degree = len(coeffs) - 1
    terms = []
    for i, coef in enumerate(coeffs):
        p = degree - i
        if abs(coef) < 1e-9:
            continue
        sign = " - " if coef < 0 else " + "
        val = abs(coef)
        val_str = f"{val:g}"

        if p == 0:
            term = f"{val_str}"
        elif p == 1:
            term = f"x" if val == 1 else f"{val_str}x"
        else:
            term = f"x^{p}" if val == 1 else f"{val_str}x^{p}"

        if not terms:
            terms.append(f"-{term}" if coef < 0 else term)
        else:
            terms.append(f"{sign}{term}")

    return "".join(terms) if terms else "0"


class HornerMultiDerivativeSolver:

    def __init__(self, n, a, c):
        self.n = n  # Bậc của đa thức
        self.a = a  # Danh sách hệ số: a[i] ứng với x^i (từ bậc 0 đến bậc n)
        self.c = c  # Điểm cần tính c
        self.log = []

    def print_and_log(self, text=""):
        print(text)
        self.log.append(text)

    def solve(self):
        self.print_and_log("=" * 70)
        self.print_and_log(
            f"TÍNH ĐẠO HÀM CÁC CẤP CỦA ĐA THỨC TẠI c = {self.c} (SƠ ĐỒ HORNER LẶP)"
        )
        self.print_and_log("=" * 70)

        poly_coeffs = [self.a[i] for i in range(self.n, -1, -1)]
        self.print_and_log(f"• Đa thức ban đầu: P(x) = {format_poly(poly_coeffs)}")
        self.print_and_log(f"• Bậc đa thức: n = {self.n}")
        self.print_and_log(f"• Điểm tính: c = {self.c}\n")

        self.print_and_log("Cơ sở lý thuyết:")
        self.print_and_log(
            "  Thực hiện chia Horner liên tiếp cho (x - c) để khai triển Taylor:"
        )
        self.print_and_log(
            "  P(x) = b_0 + b_1*(x - c) + b_2*(x - c)^2 + ... + b_n*(x - c)^n"
        )
        self.print_and_log(
            "  Từ đó đạo hàm cấp k tại c là: P^(k)(c) = k! * b_k\n"
        )

        # Ma trận lưu lại các bước Horner để in bảng
        # current_coeffs lưu các hệ số từ bậc cao xuống bậc thấp
        current_coeffs = [self.a[i] for i in range(self.n, -1, -1)]
        b = [0.0] * (self.n + 1)
        table_history = []
        table_history.append(("Hệ số ban đầu", list(current_coeffs)))

        # Horner lặp n+1 lần để tìm các hệ số b_0, b_1, ..., b_n
        for step in range(self.n + 1):
            deg = len(current_coeffs) - 1
            if deg == 0:
                b[step] = current_coeffs[0]
                break

            next_coeffs = [0.0] * deg
            prev = current_coeffs[0]
            next_coeffs[0] = prev

            for j in range(1, deg):
                prev = prev * self.c + current_coeffs[j]
                next_coeffs[j] = prev

            # Phần dư của bước chia này là b[step]
            remainder = prev * self.c + current_coeffs[deg]
            b[step] = remainder

            # Lưu lại dòng vừa tính
            row_display = list(next_coeffs) + [remainder]
            table_history.append((f"Lần chia {step + 1}", row_display))

            # Thương của bước này trở thành đa thức bị chia cho bước tiếp theo
            current_coeffs = next_coeffs

        # In bảng Horner tổng quát
        self.print_and_log("--- BẢNG HORNER LIÊN TIẾP ---")
        col_width = 12
        header = f"{'Bước':<16} | " + " | ".join(
            f"{'x^' + str(self.n - i):>{col_width}}" for i in range(self.n + 1)
        )
        self.print_and_log("-" * len(header))
        self.print_and_log(header)
        self.print_and_log("-" * len(header))

        for label, row in table_history:
            # Canh lề các số
            row_str = " | ".join(f"{val:>{col_width}.4f}" for val in row)
            self.print_and_log(f"{label:<16} | {row_str}")
        self.print_and_log("-" * len(header))

        # Tính toán giá trị đạo hàm các cấp
        self.print_and_log("\n--- KẾT QUẢ ĐẠO HÀM CÁC CẤP ---")
        self.print_and_log(
            f"• Cấp 0 (Giá trị P({self.c})):  b_0 = {b[0]:.4f} => P({self.c}) = {b[0]:.4f}"
        )

        results = []
        for k in range(1, self.n + 1):
            fact = math.factorial(k)
            deriv_val = fact * b[k]
            results.append((k, b[k], fact, deriv_val))
            self.print_and_log(
                f"• Đạo hàm cấp {k}: P^({k})({self.c}) = {k}! * b_{k} "
                f"= {fact} * ({b[k]:.4f}) = {deriv_val:.4f}"
            )

        if self.n == 0:
            self.print_and_log("Đa thức bậc 0 nên không có đạo hàm từ cấp 1 trở lên.")

        self.print_and_log("\n" + "=" * 70)
        self.print_and_log(f"TỔNG KẾT CÁC ĐẠO HÀM TẠI c = {self.c}:")
        for k, _, _, val in results:
            prime_marks = "'" * k if k <= 3 else f"^({k})"
            self.print_and_log(f"   P{prime_marks}({self.c}) = {val:.6f}")
        self.print_and_log("=" * 70 + "\n")


def main():
    try:
        with open(INP_PATH, "r", encoding="utf-8") as f:
            tokens = f.read().split()
    except FileNotFoundError:
        print(f"Lỗi: Không tìm thấy file {INP_PATH}!")
        return

    if not tokens:
        print(f"Lỗi: File {INP_PATH} rỗng!")
        return

    try:
        it = iter(tokens)
        n = int(next(it))  # Bậc n

        # Đọc n+1 hệ số từ bậc cao xuống thấp (a_n đến a_0)
        a = [0.0] * (n + 1)
        for i in range(n, -1, -1):
            a[i] = float(next(it))

        # Đọc giá trị x (hoặc c)
        c = float(next(it))
    except (StopIteration, ValueError) as e:
        print(f"Lỗi: Dữ liệu trong file {INP_PATH} không đúng định dạng!")
        return

    # Khởi tạo bộ giải và thực thi
    solver = HornerMultiDerivativeSolver(n, a, c)
    solver.solve()

    # Ghi toàn bộ kết quả vào outdaoham.txt
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(solver.log))

    print(f">>> Kết quả đã được ghi vào file: {OUT_PATH}")


if __name__ == "__main__":
    main()