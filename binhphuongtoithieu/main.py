import numpy as np
import matplotlib.pyplot as plt
import contextlib
import sys
from pathlib import Path


class Tee:
    """Ghi cùng một nội dung ra màn hình và file kết quả."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, noi_dung):
        for stream in self.streams:
            stream.write(noi_dung)
            stream.flush()
        return len(noi_dung)

    def flush(self):
        for stream in self.streams:
            stream.flush()


# ============================================================
# ĐỌC DỮ LIỆU
# ============================================================

def doc_input(ten_file):
    """
    Cho phép đọc dữ liệu theo hai dạng:

    Dạng ngang:
        Dòng 1: x1 x2 x3 ...
        Dòng 2: y1 y2 y3 ...

    Dạng dọc:
        x1 y1
        x2 y2
        x3 y3
    """

    print("=" * 60)
    print("KIEU DU LIEU DAU VAO")
    print("1. Dang ngang")
    print("2. Dang doc")

    while True:
        try:
            chon = int(input("Chon kieu dau vao: "))
            if chon in [1, 2]:
                break
            print("Vui long chon 1 hoac 2.")
        except ValueError:
            print("Gia tri khong hop le.")

    try:
        if chon == 1:
            with open(ten_file, "r", encoding="utf-8") as inp:
                dong_x = inp.readline().strip()
                dong_y = inp.readline().strip()

            if dong_x == "":
                raise ValueError("Khong co du lieu x.")

            x = np.array(dong_x.split(), dtype=float)

            # Trường hợp chỉ có x, tự sinh y = f(x)
            if dong_y == "":
                print("Chua co du lieu y.")
                print("Chuong trinh se sinh y theo ham duoc chon sau.")
                y = None
            else:
                y = np.array(dong_y.split(), dtype=float)

        else:
            x = []
            y = []

            with open(ten_file, "r", encoding="utf-8") as inp:
                for so_dong, dong in enumerate(inp, start=1):
                    dong = dong.strip()

                    if dong == "":
                        continue

                    # split() xử lý được cả dấu cách và tab
                    phan_tu = dong.split()

                    if len(phan_tu) != 2:
                        raise ValueError(
                            f"Dong {so_dong} phai co dung 2 gia tri x va y."
                        )

                    x.append(float(phan_tu[0]))
                    y.append(float(phan_tu[1]))

            x = np.array(x, dtype=float)
            y = np.array(y, dtype=float)

        return x, y

    except FileNotFoundError:
        print(f"Khong tim thay file: {ten_file}")
        raise
    except ValueError as loi:
        print("Loi du lieu:", loi)
        raise


# ============================================================
# KIỂM TRA DỮ LIỆU
# ============================================================

def kiem_tra_input(x, y):
    """
    Kiểm tra:
    - x và y phải có cùng số phần tử.
    - x không được trùng nhau.
    - x, y phải là số hữu hạn.
    """

    if y is None:
        print("Chua co du lieu y de kiem tra.")
        return 0

    if x.ndim != 1 or y.ndim != 1:
        print("Du lieu x va y phai la vector mot chieu.")
        return 0

    if x.shape[0] != y.shape[0]:
        print("Kich thuoc x va y khong hop le.")
        return 0

    if x.shape[0] < 2:
        print("Can it nhat 2 diem du lieu.")
        return 0

    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        print("Du lieu chua gia tri khong hop le.")
        return 0

    gia_tri_x, so_lan = np.unique(x, return_counts=True)

    if np.any(so_lan > 1):
        gia_tri_trung = gia_tri_x[so_lan > 1]
        print("Du lieu x bi trung:", gia_tri_trung)
        return 0

    print("Input hop le.")
    return 1


# ============================================================
# CÁC HÀM CƠ SỞ PHI
# ============================================================

def phi_x(x):
    return x


def phi_x2(x):
    return x ** 2


def phi_x3(x):
    return x ** 3


def phi_x4(x):
    return x ** 4


def phi_sin_x(x):
    return np.sin(x)


def phi_cos_x(x):
    return np.cos(x)


def phi_sin_2x(x):
    return np.sin(2 * x)


def phi_cos_2x(x):
    return np.cos(2 * x)


def phi_exp_x(x):
    return np.exp(x)


def phi_exp_am_x(x):
    return np.exp(-x)


def chon_co_so():
    """
    Chọn các hàm phi dùng để xây dựng mô hình.

    Ma trận M luôn có dạng:

        M = [1, phi1(x), phi2(x), ...]

    Phần tử đầu tiên 1 tương ứng với hệ số tự do.
    """

    danh_sach_phi = [
        ("x", phi_x),
        ("x^2", phi_x2),
        ("x^3", phi_x3),
        ("x^4", phi_x4),
        ("sin(x)", phi_sin_x),
        ("cos(x)", phi_cos_x),
        ("sin(2x)", phi_sin_2x),
        ("cos(2x)", phi_cos_2x),
        ("exp(x)", phi_exp_x),
        ("exp(-x)", phi_exp_am_x),
    ]

    print("\n" + "=" * 60)
    print("CHON CAC HAM CO SO PHI")
    print("Co the chon nhieu ham, ngan cach boi dau phay.")
    print("Vi du: 1,5,6")
    print("-" * 60)

    for i, (ten, _) in enumerate(danh_sach_phi, start=1):
        print(f"{i}. {ten}")

    while True:
        try:
            lua_chon = input("Nhap lua chon: ")
            chi_so = [int(i.strip()) for i in lua_chon.split(",")]

            if len(chi_so) == 0:
                print("Phai chon it nhat mot ham.")
                continue

            if any(i < 1 or i > len(danh_sach_phi) for i in chi_so):
                print("Chi so khong hop le.")
                continue

            if len(set(chi_so)) != len(chi_so):
                print("Khong duoc chon trung ham.")
                continue

            break

        except ValueError:
            print("Nhap theo dang: 1,5,6")

    ten_phi = [danh_sach_phi[i - 1][0] for i in chi_so]
    phi = [danh_sach_phi[i - 1][1] for i in chi_so]

    print("Da chon:", ", ".join(ten_phi))
    print("Phuong trinh cac ham co so:")
    print(" + ".join(
        f"phi{i}(x) = {ten}"
        for i, ten in enumerate(ten_phi, start=1)
    ))
    return phi, ten_phi


def in_phuong_trinh_tong_quat(kieu, ten_phi):
    """Hiển thị phương trình mô hình sau khi đã chọn các hàm cơ sở."""

    if kieu == 1:
        cac_hang_tu = " + ".join(
            f"b{i}*{ten}"
            for i, ten in enumerate(ten_phi, start=1)
        )
        print(f"Phuong trinh da chon: y = A * exp({cac_hang_tu})")
    elif kieu == 3:
        cac_hang_tu = " + ".join(
            f"a{i}*{ten}"
            for i, ten in enumerate(ten_phi, start=1)
        )
        print(f"Phuong trinh da chon: y = ln(a0 + {cac_hang_tu})")
    else:
        cac_hang_tu = " + ".join(
            f"a{i}*{ten}"
            for i, ten in enumerate(ten_phi, start=1)
        )
        print(f"Phuong trinh da chon: y = a0 + {cac_hang_tu}")


# ============================================================
# XÂY DỰNG MA TRẬN M
# ============================================================

def M(x, phi):
    """
    Xây dựng ma trận:

        M(x) = [1, phi1(x), phi2(x), ...]

    Nếu x có n phần tử và có k hàm cơ sở thì M có kích thước
    n x (k + 1).
    """

    x = np.asarray(x, dtype=float)

    cac_cot = [np.ones_like(x)]

    for ham_phi in phi:
        cac_cot.append(np.asarray(ham_phi(x), dtype=float))

    return np.column_stack(cac_cot)


# ============================================================
# BÌNH PHƯƠNG TỐI THIỂU
# ============================================================

def bpcu_tieu(x, y, phi):
    """
    Tính hệ số a theo phương pháp bình phương tối thiểu:

        M.T @ M @ a = M.T @ y

    Sử dụng np.linalg.solve.
    Nếu ma trận suy biến, dùng np.linalg.lstsq.
    """

    ma_tran_M = M(x, phi)
    MT_M = ma_tran_M.T @ ma_tran_M
    MT_y = ma_tran_M.T @ y

    print("=" * 60)
    print("MA TRAN M")
    print(ma_tran_M)

    print("=" * 60)
    print("MA TRAN M.T @ M")
    print(MT_M)

    print("=" * 60)
    print("MA TRAN M.T @ y")
    print(MT_y)

    try:
        he_so = np.linalg.solve(MT_M, MT_y)
    except np.linalg.LinAlgError:
        print("=" * 60)
        print("Ma tran M.T @ M suy bien.")
        print("Chuyen sang dung np.linalg.lstsq().")
        he_so = np.linalg.lstsq(ma_tran_M, y, rcond=None)[0]

    print("=" * 60)
    print("VECTOR HE SO")
    print(he_so)

    return he_so


def gia_tri_du_doan(x, phi, he_so):
    """
    Tính:

        y = M(x) @ he_so
    """

    return M(np.asarray(x), phi) @ he_so


def tinh_rmse(y_thuc, y_du_doan):
    return np.sqrt(np.mean((y_thuc - y_du_doan) ** 2))


# ============================================================
# CÁC MÔ HÌNH HÀM
# ============================================================

def ham_mu(x, y, phi):
    """
    Mô hình:

        y = sgn * exp(a0 + a1*phi1(x) + ...)

    Sau khi lấy log:

        ln(|y|) = a0 + a1*phi1(x) + ...

    Hệ số ban đầu của mô hình là:

        A = sgn * exp(a0)
    """

    if np.all(y > 0):
        dau = 1
    elif np.all(y < 0):
        dau = -1
    else:
        raise ValueError(
            "Ham mu yeu cau tat ca y phai cung dau va khac 0."
        )

    he_so_log = bpcu_tieu(x, np.log(np.abs(y)), phi)

    return he_so_log, dau


def ham_luy_thua(x, y):
    """
    Mô hình:

        y = A * x^b

    Lấy log:

        ln(|y|) = ln(|A|) + b*ln(x)

    Yêu cầu x > 0 và y cùng dấu.
    """

    if not np.all(x > 0):
        raise ValueError("Ham luy thua yeu cau x > 0.")

    if np.all(y > 0):
        dau = 1
    elif np.all(y < 0):
        dau = -1
    else:
        raise ValueError(
            "Ham luy thua yeu cau tat ca y phai cung dau va khac 0."
        )

    X_log = np.log(x)
    Y_log = np.log(np.abs(y))

    # Mô hình ln(|y|) = c0 + c1*ln(x)
    ma_tran = np.column_stack((np.ones_like(X_log), X_log))
    print("=" * 60)
    print("BUOC 1: CHUYEN DOI MO HINH LUY THUA")
    print("ln(|y|) = ln(|A|) + b*ln(x)")
    print("Ma tran du lieu sau khi lay log:")
    print(ma_tran)
    print("Vector ln(|y|):")
    print(Y_log)

    he_so_log = np.linalg.lstsq(ma_tran, Y_log, rcond=None)[0]
    print("BUOC 2: GIAI HE BINH PHUONG TOI THIEU")
    print("He so [ln(|A|), b]:")
    print(he_so_log)

    return he_so_log, dau


def ham_ln(x, y, phi):
    """
    Mô hình:

        y = ln(a0 + a1*phi1(x) + ...)

    Suy ra:

        exp(y) = a0 + a1*phi1(x) + ...

    """

    return bpcu_tieu(x, np.exp(y), phi)


# ============================================================
# VẼ ĐỒ THỊ
# ============================================================

def ve_do_thi(x, y, ham_du_doan, tieu_de):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x1 = np.linspace(np.min(x), np.max(x), 1000)
    y1 = ham_du_doan(x1)

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, color="red", label="Du lieu ban dau")
    plt.plot(x1, y1, color="blue", label="Duong xap xi")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(tieu_de)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================
# IN THÔNG TIN MÔ HÌNH
# ============================================================

def in_mo_hinh_tuyen_tinh(he_so, ten_phi):
    print("\nMo hinh xap xi:")

    danh_sach = [f"{he_so[0]:.10g}"]

    for i, ten in enumerate(ten_phi, start=1):
        danh_sach.append(f"({he_so[i]:.10g})*{ten}")

    print("y =", " + ".join(danh_sach))


# ============================================================
# DỰ ĐOÁN GIÁ TRỊ
# ============================================================

def du_doan_tuong_tac(ham_du_doan):
    print("\nNhap q de ket thuc.")

    while True:
        gia_tri = input("Nhap x can tinh: ").strip()

        if gia_tri.lower() == "q":
            break

        try:
            x0 = float(gia_tri)
            y0 = ham_du_doan(np.array([x0]))[0]
            print(f"y({x0}) = {y0}")
        except ValueError as loi:
            print("Khong the tinh:", loi)


# ============================================================
# CHƯƠNG TRÌNH CHÍNH
# ============================================================

def main():
    thu_muc_chuong_trinh = Path(__file__).resolve().parent
    ten_output = thu_muc_chuong_trinh / "output.txt"

    ten_file = input(
        "Nhap ten file du lieu [input1.txt]: "
    ).strip()

    if ten_file == "":
        ten_file = str(thu_muc_chuong_trinh / "input1.txt")
    elif not Path(ten_file).is_absolute():
        ten_file = str(thu_muc_chuong_trinh / ten_file)

    x, y = doc_input(ten_file)

    # Nếu y chưa có thì yêu cầu nhập hàm sinh y
    if y is None:
        print(
            "\nChuong trinh hien tai yeu cau file phai co ca x va y."
        )
        return

    if kiem_tra_input(x, y) == 0:
        return

    print("\n" + "=" * 60)
    print("CHON KIEU HAM SO")
    print("1. y = A * exp(b1*phi1(x) + b2*phi2(x) + ...)")
    print("2. y = A * x^b")
    print("3. y = ln(a0 + a1*phi1(x) + ...)")
    print("4. y = a0 + a1*phi1(x) + ...")

    while True:
        try:
            kieu = int(input("Chon kieu ham so: "))
            if kieu in [1, 2, 3, 4]:
                break
            print("Vui long chon tu 1 den 4.")
        except ValueError:
            print("Gia tri khong hop le.")

    try:
        with open(ten_output, "w", encoding="utf-8") as output_file:
            output_file.write("KET QUA TINH BINH PHUONG TOI THIEU\n")
            output_file.write(f"File du lieu: {ten_file}\n")
            output_file.write(f"Kieu ham so: {kieu}\n")
            output_file.write(f"x = {x}\n")
            output_file.write(f"y = {y}\n\n")

            # Giữ nguyên việc hiển thị trên màn hình, đồng thời lưu toàn bộ
            # các bước tính toán vào output.txt.
            with contextlib.redirect_stdout(Tee(sys.stdout, output_file)):
                if kieu == 1:
                    phi, ten_phi = chon_co_so()
                    in_phuong_trinh_tong_quat(kieu, ten_phi)
                    he_so_log, dau = ham_mu(x, y, phi)

                    A = dau * np.exp(he_so_log[0])
                    he_so_mu = he_so_log[1:]

                    def ham_du_doan(x_moi):
                        return dau * np.exp(
                            gia_tri_du_doan(x_moi, phi, he_so_log)
                        )

                    y_hat = ham_du_doan(x)
                    rmse = tinh_rmse(y, y_hat)

                    print("\n" + "=" * 60)
                    print("KET QUA HAM MU")
                    print(f"A = {A}")
                    print(f"Cac he so b = {he_so_mu}")
                    print(f"RMSE = {rmse}")

                    ve_do_thi(
                        x,
                        y,
                        ham_du_doan,
                        "Xap xi binh phuong toi thieu - Ham mu"
                    )

                    du_doan_tuong_tac(ham_du_doan)

                elif kieu == 2:
                    he_so_log, dau = ham_luy_thua(x, y)

                    A = dau * np.exp(he_so_log[0])
                    b = he_so_log[1]

                    def ham_du_doan(x_moi):
                        x_moi = np.asarray(x_moi, dtype=float)

                        if np.any(x_moi <= 0):
                            raise ValueError("Ham luy thua yeu cau x > 0.")

                        return A * x_moi ** b

                    y_hat = ham_du_doan(x)
                    rmse = tinh_rmse(y, y_hat)

                    print("\n" + "=" * 60)
                    print("KET QUA HAM LUY THUA")
                    print(f"A = {A}")
                    print(f"b = {b}")
                    print(f"RMSE = {rmse}")
                    print(f"Mo hinh: y = {A} * x^({b})")

                    ve_do_thi(
                        x,
                        y,
                        ham_du_doan,
                        "Xap xi binh phuong toi thieu - Ham luy thua"
                    )

                    du_doan_tuong_tac(ham_du_doan)

                elif kieu == 3:
                    phi, ten_phi = chon_co_so()
                    in_phuong_trinh_tong_quat(kieu, ten_phi)
                    he_so = ham_ln(x, y, phi)

                    def ham_du_doan(x_moi):
                        gia_tri_ben_trong = gia_tri_du_doan(
                            x_moi,
                            phi,
                            he_so
                        )

                        if np.any(gia_tri_ben_trong <= 0):
                            raise ValueError(
                                "Bieu thuc ben trong ln phai lon hon 0."
                            )

                        return np.log(gia_tri_ben_trong)

                    y_hat = ham_du_doan(x)
                    rmse = tinh_rmse(y, y_hat)

                    print("\n" + "=" * 60)
                    print("KET QUA HAM LOGARIT")
                    print(f"Cac he so = {he_so}")
                    print(f"RMSE = {rmse}")

                    in_mo_hinh_tuyen_tinh(he_so, ten_phi)

                    ve_do_thi(
                        x,
                        y,
                        ham_du_doan,
                        "Xap xi binh phuong toi thieu - Ham logarit"
                    )

                    du_doan_tuong_tac(ham_du_doan)

                else:
                    phi, ten_phi = chon_co_so()
                    in_phuong_trinh_tong_quat(kieu, ten_phi)
                    he_so = bpcu_tieu(x, y, phi)

                    def ham_du_doan(x_moi):
                        return gia_tri_du_doan(x_moi, phi, he_so)

                    y_hat = ham_du_doan(x)
                    rmse = tinh_rmse(y, y_hat)

                    print("\n" + "=" * 60)
                    print("KET QUA HAM TUYEN TINH TO HOP")
                    print(f"Cac he so = {he_so}")
                    print(f"RMSE = {rmse}")

                    in_mo_hinh_tuyen_tinh(he_so, ten_phi)

                    ve_do_thi(
                        x,
                        y,
                        ham_du_doan,
                        "Xap xi binh phuong toi thieu - Ham tuyen tinh"
                    )

                    du_doan_tuong_tac(ham_du_doan)

    except ValueError as loi:
        print("\nKhong the thuc hien:", loi)


if __name__ == "__main__":
    main()