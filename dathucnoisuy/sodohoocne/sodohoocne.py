import sys


def format_poly(coeffs):
    """Hàm hỗ trợ in đa thức đẹp mắt từ danh sách hệ số [a_m, a_{m-1}, ..., a_0]"""
    degree = len(coeffs) - 1
    terms = []
    for i, coef in enumerate(coeffs):
        p = degree - i
        if abs(coef) < 1e-9:
            continue
        # Dấu và giá trị tuyệt đối
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


class HornerSolver:

    def __init__(self, n, a, c):
        self.n = n  # Bậc của đa thức
        self.a = a  # Danh sách hệ số: a[i] ứng với x^i
        self.c = c  # Giá trị x = c
        self.log = []

    def print_and_log(self, text=""):
        print(text)
        self.log.append(text)

    # -------------------------------------------------------------
    # CASE 1: TÍNH GIÁ TRỊ ĐA THỨC THEO SƠ ĐỒ HORNER
    # -------------------------------------------------------------
    def case1_tinh_gia_tri(self):
        self.print_and_log("=" * 60)
        self.print_and_log(
            f"CASE 1: TÍNH GIÁ TRỊ ĐA THỨC P({self.c}) BẰNG SƠ ĐỒ HORNER"
        )
        self.print_and_log("=" * 60)

        poly_coeffs = [self.a[i] for i in range(self.n, -1, -1)]
        self.print_and_log(f"Đa thức: P(x) = {format_poly(poly_coeffs)}")
        self.print_and_log(
            "Công thức truy hồi Horner: b_n = a_n;  b_k = b_{k+1} * c + a_k (với k = n-1 down to 0)\n"
        )

        b = [0.0] * (self.n + 1)
        b[self.n] = self.a[self.n]
        self.print_and_log(f"• Bước k = {self.n}:  b_{self.n} = a_{self.n} = {b[self.n]}")

        for k in range(self.n - 1, -1, -1):
            b[k] = b[k + 1] * self.c + self.a[k]
            self.print_and_log(
                f"• Bước k = {k:2d}: b_{k} = b_{k+1} * ({self.c}) + a_{k} "
                f"= {b[k+1]} * ({self.c}) + ({self.a[k]}) = {b[k]:.4f}"
            )

        self.print_and_log(f"\n=> ĐÁP ÁN: P({self.c}) = b_0 = {b[0]:.4f}\n")
        return b

    # -------------------------------------------------------------
    # CASE 2: TÍNH GIÁ TRỊ ĐẠO HÀM P'(c)
    # -------------------------------------------------------------
    def case2_tinh_dao_ham(self):
        self.print_and_log("=" * 60)
        self.print_and_log(f"CASE 2: TÍNH GIÁ TRỊ ĐẠO HÀM P'({self.c})")
        self.print_and_log("=" * 60)

        if self.n == 0:
            self.print_and_log("Đa thức bậc 0 nên đạo hàm P'(x) = 0")
            self.print_and_log(f"=> ĐÁP ÁN: P'({self.c}) = 0\n")
            return

        # Tính các hệ số đạo hàm: d_{i-1} = i * a_i
        d = [0.0] * self.n
        self.print_and_log("1. Xác định các hệ số của đa thức đạo hàm P'(x):")
        for i in range(self.n, 0, -1):
            d[i - 1] = i * self.a[i]
            self.print_and_log(
                f"   d_{i-1} = {i} * a_{i} = {i} * ({self.a[i]}) = {d[i-1]}"
            )

        deriv_coeffs = [d[i] for i in range(self.n - 1, -1, -1)]
        self.print_and_log(
            f"\n   Đa thức đạo hàm: P'(x) = {format_poly(deriv_coeffs)}"
        )

        self.print_and_log(
            f"\n2. Áp dụng sơ đồ Horner tính giá trị P'({self.c}):"
        )
        b_deriv = [0.0] * self.n
        b_deriv[self.n - 1] = d[self.n - 1]
        self.print_and_log(
            f"• Bước k = {self.n-1}: b_{self.n-1} = d_{self.n-1} = {b_deriv[self.n-1]}"
        )

        for k in range(self.n - 2, -1, -1):
            b_deriv[k] = b_deriv[k + 1] * self.c + d[k]
            self.print_and_log(
                f"• Bước k = {k:2d}: b_{k} = b_{k+1} * ({self.c}) + d_{k} "
                f"= {b_deriv[k+1]} * ({self.c}) + ({d[k]}) = {b_deriv[k]:.4f}"
            )

        self.print_and_log(
            f"\n=> ĐÁP ÁN: Đạo hàm P'({self.c}) = {b_deriv[0]:.4f}\n"
        )

    # -------------------------------------------------------------
    # CASE 3: CHIA ĐA THỨC CHO (x - c)
    # -------------------------------------------------------------
    def case3_chia_da_thuc(self):
        self.print_and_log("=" * 60)
        self.print_and_log(f"CASE 3: CHIA ĐA THỨC P(x) CHO NHỊ THỨC (x - {self.c})")
        self.print_and_log("=" * 60)
        self.print_and_log("Sơ đồ Horner: P(x) = (x - c) * Q(x) + R")

        b = [0.0] * (self.n + 1)
        b[self.n] = self.a[self.n]
        for k in range(self.n - 1, -1, -1):
            b[k] = b[k + 1] * self.c + self.a[k]

        # Bảng tóm tắt hệ số
        self.print_and_log("\nBảng tính Horner:")
        row_a = f"{'Hệ số a_k':<12} | " + " | ".join(
            f"{self.a[i]:>8.2f}" for i in range(self.n, -1, -1)
        )
        row_b = f"{'Hệ số b_k':<12} | " + " | ".join(
            f"{b[i]:>8.2f}" for i in range(self.n, -1, -1)
        )
        self.print_and_log("-" * len(row_a))
        self.print_and_log(row_a)
        self.print_and_log(row_b)
        self.print_and_log("-" * len(row_a))

        q_coeffs = [b[i] for i in range(self.n, 0, -1)]
        du = b[0]

        self.print_and_log(
            f"\n• Các hệ số của thương Q(x): {[round(x, 4) for x in q_coeffs]}"
        )
        self.print_and_log(f"• Phần dư R = {du:.4f}")
        self.print_and_log(
            f"\n=> ĐÁP ÁN: P(x) = (x - {self.c}) * [{format_poly(q_coeffs)}] + ({du:.4f})\n"
        )

    # -------------------------------------------------------------
    # CASE 4: NHÂN ĐA THỨC VỚI (x - c)
    # -------------------------------------------------------------
    def case4_nhan_da_thuc(self):
        self.print_and_log("=" * 60)
        self.print_and_log(f"CASE 4: NHÂN ĐA THỨC P(x) VỚI NHỊ THỨC (x - {self.c})")
        self.print_and_log("=" * 60)
        self.print_and_log(
            "Công thức: P(x) * (x - c) có bậc n + 1.\nCác hệ số C_k của x^k được tính như sau:"
        )

        # Đa thức tích có bậc n+1, mảng hệ số từ bậc 0 đến n+1
        res = [0.0] * (self.n + 2)

        # C_{n+1} = a_n
        res[self.n + 1] = self.a[self.n]
        self.print_and_log(
            f"• Bậc x^{self.n+1}: C_{self.n+1} = a_{self.n} = {res[self.n+1]}"
        )

        # C_k = a_{k-1} - c * a_k (với k từ n giảm về 1)
        for k in range(self.n, 0, -1):
            res[k] = self.a[k - 1] - self.c * self.a[k]
            self.print_and_log(
                f"• Bậc x^{k:2d}: C_{k} = a_{k-1} - ({self.c}) * a_{k} "
                f"= {self.a[k-1]} - ({self.c}) * ({self.a[k]}) = {res[k]:.4f}"
            )

        # C_0 = -c * a_0
        res[0] = -self.c * self.a[0]
        self.print_and_log(
            f"• Bậc x^0 : C_0 = -({self.c}) * a_0 = -({self.c}) * ({self.a[0]}) = {res[0]:.4f}"
        )

        prod_coeffs = [res[i] for i in range(self.n + 1, -1, -1)]
        self.print_and_log(
            f"\n=> ĐÁP ÁN: P(x) * (x - {self.c}) = {format_poly(prod_coeffs)}\n"
        )


def main():
    try:
        with open("inp.txt", "r") as f:
            tokens = f.read().split()
    except FileNotFoundError:
        print("Lỗi: Không tìm thấy file inp.txt!")
        return

    if not tokens:
        print("File inp.txt rỗng!")
        return

    it = iter(tokens)
    n = int(next(it))  # Bậc n

    # Đọc n+1 hệ số từ bậc cao xuống thấp (a_n đến a_0)
    a = [0.0] * (n + 1)
    for i in range(n, -1, -1):
        a[i] = float(next(it))

    # Đọc giá trị x (hoặc c)
    x = float(next(it))

    # Khởi tạo bộ giải và thực thi cả 4 case
    solver = HornerSolver(n, a, x)
    solver.case1_tinh_gia_tri()
    solver.case2_tinh_dao_ham()
    solver.case3_chia_da_thuc()
    solver.case4_nhan_da_thuc()

    # Lưu toàn bộ lời giải và kết quả vào out.txt
    with open("out.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(solver.log))

    print(">>> Đã thực hiện xong cả 4 case và ghi chi tiết vào out.txt!")


if __name__ == "__main__":
    main()