import os
import sys
import numpy as np

log_file_handle = None


def pr(text=""):
    print(text)
    if log_file_handle:
        print(text, file=log_file_handle)


def doc_input(ten_file):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, ten_file)

    if not os.path.exists(file_path):
        print(f"Lỗi: Không tìm thấy file {file_path}")
        sys.exit()

    with open(file_path, "r", encoding="utf-8") as inp:
        lines = [line.strip() for line in inp.readlines() if line.strip()]

    x = np.array(lines[0].split(), dtype=float)
    y = np.array(lines[1].split(), dtype=float)
    return x, y


def hoocne_product(x):
    a = []
    Li = [np.array([1.0])]
    a.append(np.array([1.0, 0.0]))
    for root in x:
        b = a[-1]
        c = [1.0]
        for j in range(len(b) - 1):
            c.append(b[j + 1] - b[j] * root)
        c.append(0.0)
        a.append(np.array(c))
        Li.append(np.delete(a[-1], -1))
    return Li


def hoocne_quatient(a, x0):
    y = [a[0]]
    for i in range(len(a) - 1):
        y.append(y[i] * x0 + a[i + 1])
    b = np.array(y[:-1])
    b_0 = y[-1]
    return b, b_0


# TÌM TẤT CẢ CÁC KHOẢNG ĐƠN ĐIỆU (SONG ÁNH) CỦA DỮ LIỆU
def tim_tat_ca_khoang_song_anh(x, y):
    n = len(x)
    sgn = 1 if y[0] < y[1] else -1
    danh_sach = []
    start_idx = 0

    for i in range(1, n):
        curr_sgn = 1 if y[i - 1] < y[i] else -1
        if curr_sgn != sgn:
            danh_sach.append((start_idx, i - 1, sgn))
            start_idx = i - 1
            sgn = curr_sgn
    danh_sach.append((start_idx, n - 1, sgn))
    return danh_sach


# LẶP NEWTON TIẾN
def lap_newton_tien(x, y, y0, eps):
    pr("\n" + "=" * 70)
    pr("               CHI TIẾT PHƯƠNG PHÁP LẶP NEWTON TIẾN")
    pr("=" * 70)
    n = len(y)
    h = x[1] - x[0]
    deltaY0 = y[1] - y[0]

    pr(f"• Mốc gốc x_0 = {x[0]}, y_0 = {y[0]}")
    pr(f"• Bước nhảy h = {h:.6f}")
    pr(f"• Sai phân cấp 1: Δy_0 = {deltaY0:.6f}")

    t0 = (y0 - y[0]) / deltaY0
    t1 = t0
    pr(f"• Giá trị xấp xỉ ban đầu: t^(0) = ({y0} - {y[0]}) / {deltaY0:.6f} = {t0:.6f}")

    P = [np.array([y[0] - y0]), np.array([y[0] - y0, 0.0])]
    Z = [np.array([y[0]]), np.array([y[1], y[1] - y[0]])]

    t_arr = np.array(range(n))
    Li = hoocne_product(t_arr)
    giai_thua = 1
    lan_lap = 0

    pr("\nBẢNG TIẾN TRÌNH CÁC BƯỚC LẶP NEWTON TIẾN:")
    pr(f"{'Lần lặp':<10}{'t^(m)':<15}{'Sai số |t1 - t0|':<20}{'x^(m) tương ứng':<18}")
    pr("-" * 65)
    pr(f"{lan_lap:<10}{t1:<15.6f}{'-':<20}{x[0] + t1 * h:<18.8f}")

    for i in range(2, n):
        lan_lap += 1
        giai_thua *= i
        z = [y[i]]
        for j in range(i):
            z.append(z[j] - Z[i - 1][j])
        Z.append(np.array(z))
        Fi = Z[-1][-1]

        P.append(np.concatenate((np.array([0.0]), P[i - 1])) + (Fi * Li[i]) / giai_thua)
        t_prev = t1
        val_P = hoocne_quatient(P[-1], t_prev)[1]
        t1 = -1.0 / deltaY0 * val_P
        diff = abs(t1 - t_prev)
        pr(f"{lan_lap:<10}{t1:<15.6f}{diff:<20.8e}{x[0] + t1 * h:<18.8f}")
        if diff < eps:
            break

    while abs(t1 - t0) > eps and lan_lap < 100:
        lan_lap += 1
        t_prev = t1
        val_P = hoocne_quatient(P[-1], t_prev)[1]
        t1 = -1.0 / deltaY0 * val_P
        diff = abs(t1 - t_prev)
        pr(f"{lan_lap:<10}{t1:<15.6f}{diff:<20.8e}{x[0] + t1 * h:<18.8f}")

    x_res = x[0] + t1 * h
    pr("-" * 65)
    pr(f"=> KẾT QUẢ: Nghiệm x khi y = {y0} là: x = {x_res:.8f}\n")
    return x_res


# LẶP NEWTON LÙI
def lap_newton_lui(x, y, y0, eps):
    pr("\n" + "=" * 70)
    pr("               CHI TIẾT PHƯƠNG PHÁP LẶP NEWTON LÙI")
    pr("=" * 70)
    n = len(y)
    h = x[0] - x[1]
    deltaY0 = y[1] - y[0]

    pr(f"• Mốc gốc x_0 = {x[0]}, y_0 = {y[0]}")
    pr(f"• Bước nhảy h = {h:.6f}")
    t0 = (y0 - y[0]) / deltaY0
    t1 = t0

    P = [np.array([y[0] - y0]), np.array([y[0] - y0, 0.0])]
    Z = [np.array([y[0]]), np.array([y[1], y[1] - y[0]])]

    t_arr = np.array(range(n))
    Li = hoocne_product(t_arr)
    giai_thua = 1
    lan_lap = 0

    pr("\nBẢNG TIẾN TRÌNH CÁC BƯỚC LẶP NEWTON LÙI:")
    pr(f"{'Lần lặp':<10}{'t^(m)':<15}{'Sai số |t1 - t0|':<20}{'x^(m) tương ứng':<18}")
    pr("-" * 65)
    pr(f"{lan_lap:<10}{t1:<15.6f}{'-':<20}{x[0] - t1 * h:<18.8f}")

    for i in range(2, n):
        lan_lap += 1
        giai_thua *= i
        z = [y[i]]
        for j in range(i):
            z.append(z[j] - Z[i - 1][j])
        Z.append(np.array(z))
        Fi = Z[-1][-1]

        P.append(np.concatenate((np.array([0.0]), P[i - 1])) + (Fi * Li[i]) / giai_thua)
        t_prev = t1
        t1 = -1.0 / deltaY0 * hoocne_quatient(P[-1], t_prev)[1]
        diff = abs(t1 - t_prev)
        pr(f"{lan_lap:<10}{t1:<15.6f}{diff:<20.8e}{x[0] - t1 * h:<18.8f}")
        if diff < eps:
            break

    while abs(t1 - t0) > eps and lan_lap < 100:
        lan_lap += 1
        t_prev = t1
        t1 = -1.0 / deltaY0 * hoocne_quatient(P[-1], t_prev)[1]
        diff = abs(t1 - t_prev)
        pr(f"{lan_lap:<10}{t1:<15.6f}{diff:<20.8e}{x[0] - t1 * h:<18.8f}")

    x_res = x[0] - t1 * h
    pr("-" * 65)
    pr(f"=> KẾT QUẢ: Nghiệm x khi y = {y0} là: x = {x_res:.8f}\n")
    return x_res


# PHƯƠNG PHÁP HÀM NGƯỢC
def ham_nguoc(x_bar, y_bar, y0):
    pr("\n" + "=" * 70)
    pr("               CHI TIẾT PHƯƠNG PHÁP HÀM NGƯỢC x = P(y)")
    pr("=" * 70)
    n = len(y_bar)
    Li = hoocne_product(y_bar)
    P = []
    Z = []

    pr("1. BẢNG TỶ HIỆU VỚI CÁC MỐC y_i:")
    for i in range(n):
        z = [x_bar[i]]
        for j in range(i):
            val = (z[j] - Z[i - 1][j]) / (y_bar[i] - y_bar[i - j - 1])
            z.append(val)
        Z.append(np.array(z))
        Fi = Z[i][i]
        pr(f"   • Tỷ hiệu f[y_0...y_{i}] = {Fi:.8f}")
        if i == 0:
            P.append(np.array(Fi * Li[i]))
        else:
            P.append(np.concatenate((np.array([0.0]), P[i - 1])) + Fi * Li[i])

    poly_coeffs = P[-1]
    pr(f"\n2. HỆ SỐ ĐA THỨC HÀM NGƯỢC x = P(y): {list(np.round(poly_coeffs, 8))}")

    _, x_res = hoocne_quatient(poly_coeffs, y0)
    pr(f"\n3. TÍNH x = P({y0}) BẰNG LƯỢC ĐỒ HORNER: x = {x_res:.8f}")
    pr(f"\n=> KẾT QUẢ: Nghiệm x tương ứng với y = {y0} là: x = {x_res:.8f}\n")
    return x_res


# =====================================================================
# HÀM CHÍNH (MAIN)
# =====================================================================
def main():
    global log_file_handle

    current_dir = os.path.dirname(os.path.abspath(__file__))
    out_file_path = os.path.join(current_dir, "output.txt")
    log_file_handle = open(out_file_path, "w", encoding="utf-8")

    pr("=" * 75)
    pr("                  BÁO CÁO BÀI TOÁN NỘI SUY NGƯỢC")
    pr("=" * 75)

    x, y = doc_input("input.txt")
    pr("\n1. DỮ LIỆU ĐỌC ĐƯỢC TỪ INPUT.TXT:")
    pr(f"   • Số lượng mốc: {len(x)}")
    pr(f"   • Khoảng x: [{x[0]}, {x[-1]}]")
    pr(f"   • Giá trị y nhỏ nhất: {np.min(y):.5f}, lớn nhất: {np.max(y):.5f}")

    # BƯỚC 1: TÌM VÀ IN TOÀN BỘ CÁC KHOẢNG SONG ÁNH TRƯỚC
    tat_ca_khoang = tim_tat_ca_khoang_song_anh(x, y)
    pr("\n2. DANH SÁCH TẤT CẢ CÁC KHOẢNG SONG ÁNH CỦA BẢNG DỮ LIỆU:")
    print("-" * 75)
    for idx, (s, e, sign) in enumerate(tat_ca_khoang, start=1):
        chieu = "NGHỊCH BIẾN (GIẢM)" if sign == -1 else "ĐỒNG BIẾN (TĂNG)"
        y_min = min(y[s], y[e])
        y_max = max(y[s], y[e])
        pr(f"   Khoảng {idx}: Đoạn x ∈ [{x[s]:.2f}, {x[e]:.2f}] --> y ∈ [{y_min:.5f}, {y_max:.5f}] ({chieu})")
    print("-" * 75)

    # BƯỚC 2: CHO NGƯỜI DÙNG NHẬP y0
    print("\n* Lưu ý: Giá trị y0 bạn nhập phải nằm trong khoảng y của ít nhất một khoảng song ánh ở trên.")
    while True:
        try:
            y0_str = input(">> Nhập giá trị y0 cần tìm x: ").strip()
            y0 = float(y0_str)
            # Kiểm tra xem y0 có thuộc khoảng nào không
            khoang_chua_y0 = [
                (idx, k)
                for idx, k in enumerate(tat_ca_khoang, start=1)
                if min(y[k[0]], y[k[1]]) <= y0 <= max(y[k[0]], y[k[1]])
            ]
            if not khoang_chua_y0:
                print(f"Lỗi: Giá trị y0 = {y0} không thuộc bất kỳ khoảng song ánh nào ở trên! Vui lòng nhập lại.")
            else:
                break
        except ValueError:
            print("Vui lòng nhập một số thực hợp lệ!")

    pr(f"\n>>> Giá trị y* cần tìm nghiệm x: y* = {y0}")

    # BƯỚC 3: HIỂN THỊ CÁC KHOẢNG CHỨA y0 VÀ CHO CHỌN
    print("\nCác khoảng song ánh THỎA MÃN (chứa giá trị y0 vừa nhập):")
    for idx, k in khoang_chua_y0:
        s, e, sign = k
        chieu = "GIẢM" if sign == -1 else "TĂNG"
        print(f"   [{idx}] Khoảng {idx}: x ∈ [{x[s]:.2f}, {x[e]:.2f}] ({chieu})")

    if len(khoang_chua_y0) == 1:
        chon_khoang = khoang_chua_y0[0][0]
        print(f">> Tự động chọn Khoảng {chon_khoang}")
    else:
        chon_khoang = int(input(f">> Chọn khoảng song ánh muốn tìm nghiệm ({[k[0] for k in khoang_chua_y0]}): "))

    k_sel = tat_ca_khoang[chon_khoang - 1]
    s_idx, e_idx, sign_sel = k_sel
    pr(f"\n>>> Bạn đã chọn Khoảng {chon_khoang}: x ∈ [{x[s_idx]:.2f}, {x[e_idx]:.2f}]")

    # Xác định vị trí lân cận y0
    index = s_idx
    for i in range(s_idx, e_idx):
        if min(y[i], y[i + 1]) <= y0 <= max(y[i], y[i + 1]):
            index = i
            break

    # BƯỚC 4: CHỌN PHƯƠNG PHÁP
    print("\nChọn phương pháp giải:")
    print("1. Lặp Newton")
    print("2. Phương pháp hàm ngược")
    print("3. Thực hiện CẢ HAI phương pháp để so sánh")
    choice = input("Nhập lựa chọn (1/2/3, mặc định 3): ").strip()
    if choice not in ["1", "2", "3"]:
        choice = "3"

    # THỰC THI LẶP NEWTON
    if choice in ["1", "3"]:
        n_total = len(x)
        eps = 1e-6
        if index < n_total / 2:
            pr("\n[Phân tích vị trí]: Điểm nằm ở nửa đầu bảng mốc -> Áp dụng LẶP NEWTON TIẾN")
            xbar = x[index:n_total]
            ybar = y[index:n_total]
            lap_newton_tien(xbar, ybar, y0, eps)
        else:
            pr("\n[Phân tích vị trí]: Điểm nằm ở nửa sau bảng mốc -> Áp dụng LẶP NEWTON LÙI")
            xbar = np.flip(x[0 : index + 2])
            ybar = np.flip(y[0 : index + 2])
            lap_newton_lui(xbar, ybar, y0, eps)

    # THỰC THI HÀM NGƯỢC
    if choice in ["2", "3"]:
        max_moc = e_idx - s_idx + 1
        num_moc_str = input(f"\nNhập số mốc dùng cho hàm ngược (2 -> {max_moc}, mặc định {max_moc}): ").strip()
        num_moc = int(num_moc_str) if (num_moc_str.isdigit() and 2 <= int(num_moc_str) <= max_moc) else max_moc

        half = num_moc // 2
        start_m = max(s_idx, index - half)
        end_m = start_m + num_moc
        if end_m > e_idx + 1:
            end_m = e_idx + 1
            start_m = end_m - num_moc

        xbar = x[start_m:end_m]
        ybar = y[start_m:end_m]
        pr(f"\nCác mốc trích xuất cho hàm ngược (tổng cộng {len(xbar)} mốc):")
        for xi, yi in zip(xbar, ybar):
            pr(f"   (x = {xi:.4f}, y = {yi:.5f})")

        ham_nguoc(xbar, ybar, y0)

    log_file_handle.close()
    print(f"\n>>> HOÀN TẤT! Đã lưu toàn bộ các bước làm vào: {out_file_path}")


if __name__ == "__main__":
    main()