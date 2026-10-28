import os
import sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(CURRENT_DIR, "out.txt")


def parse_text_tokens(raw_text):
    """Trích xuất các số thực từ chuỗi văn bản (dấu cách, tab, phẩy)."""
    # Thay dấu phẩy hoặc chấm phẩy bằng dấu cách nếu có
    cleaned_text = raw_text.replace(",", " ").replace(";", " ")
    tokens = cleaned_text.split()
    nums = []
    for t in tokens:
        try:
            nums.append(float(t))
        except ValueError:
            continue  # Bỏ qua chữ (tiêu đề cột x, y nếu có)

    points = []
    # Nếu số lượng phần tử lẻ và số đầu tiên là số lượng điểm N
    if len(nums) % 2 == 1 and len(nums) > 1:
        nums = nums[1:]  # Bỏ qua số N ở đầu

    it = iter(nums)
    while True:
        try:
            x = next(it)
            y = next(it)
            points.append((x, y))
        except StopIteration:
            break
    return points


def doc_du_lieu():
    """Hỗ trợ đọc thông minh mọi loại file: Excel thật, CSV, Text."""
    files = os.listdir(CURRENT_DIR)
    candidates = [
        f
        for f in files
        if f.endswith((".xlsx", ".xls", ".txt", ".csv", ".inp", ".dat"))
        and f != "out.txt"
        and not f.startswith("~$")
    ]

    print("=" * 60)
    print("CHỌN FILE DỮ LIỆU ĐẦU VÀO:")
    if not candidates:
        file_name = input("Nhập tên file dữ liệu: ").strip()
    else:
        print("Các file tìm thấy trong thư mục:")
        for idx, f in enumerate(candidates, 1):
            print(f"  {idx}. {f}")
        choice = input(
            f"Chọn số thứ tự [1-{len(candidates)}] (hoặc nhập tên file): "
        ).strip()
        if choice.isdigit() and 1 <= int(choice) <= len(candidates):
            file_name = candidates[int(choice) - 1]
        else:
            file_name = choice

    file_path = os.path.join(CURRENT_DIR, file_name)
    if not os.path.exists(file_path):
        print(f"Lỗi: Không tìm thấy file {file_path}!")
        sys.exit(1)

    points = []

    # 1. Thử đọc dạng Excel (chỉ định engine='openpyxl')
    doc_thanh_cong = False
    if file_name.endswith((".xlsx", ".xls")):
        try:
            import pandas as pd

            df = pd.read_excel(file_path, engine="openpyxl")
            df_clean = df.iloc[:, [0, 1]].dropna()
            for _, row in df_clean.iterrows():
                try:
                    points.append((float(row.iloc[0]), float(row.iloc[1])))
                except ValueError:
                    continue
            if len(points) > 0:
                doc_thanh_cong = True
        except Exception:
            # Nếu đọc theo dạng Excel nhị phân thất bại, chuyển sang đọc text
            doc_thanh_cong = False

    # 2. Nếu không phải Excel hoặc file .xlsx là file text đổi đuôi
    if not doc_thanh_cong:
        for enc in ["utf-8", "utf-8-sig", "latin1", "cp1252"]:
            try:
                with open(file_path, "r", encoding=enc, errors="ignore") as f:
                    content = f.read()
                points = parse_text_tokens(content)
                if len(points) > 0:
                    doc_thanh_cong = True
                    break
            except Exception:
                continue

    if not points:
        print(
            f"Lỗi: Không thể trích xuất được dữ liệu (x, y) hợp lệ từ file {file_name}!"
        )
        print("Vui lòng kiểm tra lại nội dung file.")
        sys.exit(1)

    # Sắp xếp các điểm theo thứ tự x tăng dần
    points.sort(key=lambda p: p[0])
    return file_name, points


def trich_xuat(points, x_bar, k):
    """Trích xuất k mốc tối ưu nhất quanh x_bar (chuẩn Stirling đối xứng tâm)."""
    total = len(points)
    if k > total:
        print(f"Lỗi: Bảng chỉ có {total} điểm, không thể trích xuất {k} điểm!")
        sys.exit(1)

    # Tìm mốc x_i gần x_bar nhất
    idx_center = min(range(total), key=lambda i: abs(points[i][0] - x_bar))
    x_0 = points[idx_center][0]

    m = k // 2
    if k % 2 == 1:
        start_idx = idx_center - m
        end_idx = idx_center + m
    else:
        start_idx = idx_center - m
        end_idx = idx_center + m - 1

    # Xử lý chạm mép đầu hoặc mép cuối của bảng
    if start_idx < 0:
        start_idx = 0
        end_idx = k - 1
    elif end_idx >= total:
        end_idx = total - 1
        start_idx = total - k

    selected = points[start_idx : end_idx + 1]
    return idx_center, x_0, selected


def main():
    file_name, points = doc_du_lieu()
    print(f"\n=> Đã đọc thành công {len(points)} điểm từ file: {file_name}")

    # Nhập x_bar và k
    while True:
        try:
            x_bar = float(input("\nNhập giá trị x ngang (ví dụ 7.83): "))
            break
        except ValueError:
            print("Vui lòng nhập số thực hợp lệ!")

    while True:
        try:
            k = int(
                input(
                    f"Nhập số lượng mốc cần trích xuất (tối đa {len(points)}, ví dụ 9): "
                )
            )
            if 1 <= k <= len(points):
                break
            print(f"Số mốc phải từ 1 đến {len(points)}!")
        except ValueError:
            print("Vui lòng nhập số nguyên hợp lệ!")

    idx_center, x_0, selected = trich_xuat(points, x_bar, k)

    # In kết quả
    lines_output = []
    lines_output.append("=" * 65)
    lines_output.append(f"KẾT QUẢ TRÍCH XUẤT {k} MỐC NỘI SUY TỐI ƯU")
    lines_output.append("=" * 65)
    lines_output.append(f"• File nguồn: {file_name}")
    lines_output.append(f"• Điểm x ngang: {x_bar}")
    lines_output.append(
        f"• Mốc tâm gần nhất: x_0 = {x_0} (|x_0 - x ngang| = {abs(x_0 - x_bar):.6f})\n"
    )

    header = f"{'Chỉ số':<10}{'x':<15}{'y':<15}{'Ghi chú':<15}"
    lines_output.append(header)
    lines_output.append("-" * 55)

    m = k // 2 if k % 2 == 1 else 0
    for idx, (x, y) in enumerate(selected):
        if k % 2 == 1:
            idx_str = f"i = {idx - m}"
            note = "<-- Mốc tâm x_0" if (idx - m == 0) else ""
        else:
            idx_str = f"i = {idx}"
            note = ""
        lines_output.append(f"{idx_str:<10}{x:<15.6f}{y:<15.6f}{note:<15}")

    lines_output.append("-" * 55)

    for line in lines_output:
        print(line)

    # Ghi dữ liệu ra out.txt
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        # Dòng 1: Số lượng mốc
        f.write(f"{len(selected)}\n")
        # Các dòng tiếp theo: x y
        for x, y in selected:
            f.write(f"{x} {y}\n")
        # Dòng cuối: Điểm x ngang
        f.write(f"{x_bar}\n")

    print(f"\n>>> Đã xuất dữ liệu chuẩn vào file: {OUT_PATH}")


if __name__ == "__main__":
    main()