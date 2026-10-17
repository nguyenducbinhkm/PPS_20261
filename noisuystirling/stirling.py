import numpy as np
from termcolor import colored
from math import prod
np.set_printoptions(precision=16)
m = open('interpolation.txt','r').readlines()
xs =[]
ys = []
for i in m[:-2]:
    k=i.strip().split('\t')
    xs.append(float(k[0]))
    ys.append(float(k[1]))
k = m[-1].strip().split(' ')
start = float(k[0])
def vietaFormula(coeff:list,x0:float) -> list:
    coeff = [0]+coeff
    n = len(coeff)
    for i in range(n-1):
        coeff[i]-=x0*coeff[i+1]
    return coeff # 3(n-1) mul
def plus_poly(high: list,low: list) -> list:
    n = len(low)
    for i in range(n):
        high[i]+=low[i]
    return high #len(high) cal
def mul_list(kirito: list,ratio: float) -> list:
    return [i*ratio for i in kirito] # len(kirito) cal


def deri(ys:list):
    n = len(ys)
    return [ys[i+1]-ys[i] for i in range(n-1)]    # len(ys)*2 cal
def nodes(xs:list,c:float):
    p = len(xs)
    h = xs[1]-xs[0]
    j = 0
    if abs(c-xs[0]) < abs(c-xs[-1]):
        i=0
        while abs(xs[i]-c) > h/2:
            i+=1
        if xs[i]==c:
            exit(colored(f'Giá trị tại {c} đã biết là: {ys[i]}','yellow'))
        while xs[i+j]-xs[i-j]<=1 and 0<=i+j<=p-1 and 0<=i-j<=p-1:
            j+=1
        j-=1
        if j<=1:
            if (i-2)*(i+2)>=0:
                return range(i-2,i+3)
            else:
                exit(colored('Có quá ít điểm dữ liệu để thực hiện nội suy Stirling'.upper(),'red'))
        elif 2<=j<=4:
            return range(i-j,i+j+1)
        else:
            return range(i-4,i+5)
    else:
        i=-1
        while abs(xs[i]-c) > h/2:
            i-=1
        if xs[i]==c:
            exit(colored(f'Giá trị tại {c} đã biết là: {ys[i]}','yellow'))
        while xs[i+j]-xs[i-j]<=1 and -p<=i+j<=-1 and -p<=i-j<=-1:
            j+=1
        j-=1
        if j<=1:
            if (i-2)*(i+2)>0:
                return range(i-2,i+3) # min 5 nodes
            else:
                exit(colored('Có quá ít điểm dữ liệu để thực hiện nội suy Stirling'.upper(),'red'))
        elif 2<=j<=4:
            return range(i-j,i+j+1)
        else:
            return range(i-4,i+5) # max 9 nodes
def Stirling_overhaul(xs:list,ys:list,x:float,output):
    ys = [ys[i] for i in nodes(xs,x)] ### phải cắt ys trước nếu không sẽ ảnh hưởng algorithm
    xs = [xs[i] for i in nodes(xs,x)]
    print(f'NỘI SUY STIRLING TẠI x = {x}', file=output)
    print('=' * 60, file=output)
    print('Các mốc được chọn và giá trị tương ứng:', file=output)
    for m,i in enumerate(xs):
        print(f'x = {xs[m]}, y = {ys[m]}', file=output) # in ra các điểm đã lấy
    h = round(xs[1] - xs[0], 10)
    n = int((len(xs) - 1) / 2)
    x0 = xs[n]
    y0 = ys[n]
    t = (x - x0) / h
    print(file=output)
    print(f'Bước h = {h}', file=output)
    print(f'Mốc trung tâm x0 = {x0}', file=output)
    print(f'Giá trị trung tâm y0 = {y0}', file=output)
    print(f't = (x - x0) / h = ({x} - {x0}) / {h} = {t}', file=output)
    print(file=output)
    print('Các sai phân được sử dụng:', file=output)
    coeff = [1]
    t_coeffs = [y0]
    even_coeff = []
    odd_coeff = []
    fact = 1
    for i in range(n):
      ys = deri(ys[:])

      fact *= 2 * i + 1
      a_i = (ys[n - i] + ys[n - i - 1]) / (2 * fact)
      odd_coeff = plus_poly(mul_list(coeff, a_i), odd_coeff)
      print(f'a{i + 1} = {a_i}', file=output)
      ys = deri(ys[:])

      fact *= 2 * i + 2
      b_i = ys[n - i - 1] / fact
      even_coeff = plus_poly(mul_list(coeff, b_i), even_coeff)
      print(f'b{i + 1} = {b_i}', file=output)
      coeff = vietaFormula(coeff, (i + 1) ** 2)

    odd_coeff = np.array(odd_coeff)
    even_coeff = np.array(even_coeff)
    t_sq = t**2
    q = 1
    t_vec = []
    for i in range(1, n + 1):
      q *= t_sq
      t_vec.append([q])
    t_vec = np.array(t_vec)

    result = y0 + even_coeff @ t_vec + (odd_coeff @ t_vec) / t
    print(file=output)
    print(f'Giá trị nội suy Stirling tại x = {x}:', file=output)
    print(result, file=output)

    for j in range(n):
      t_coeffs.append(float(odd_coeff[j]))
      t_coeffs.append(float(even_coeff[j]))

    print(file=output)
    print(
        f'Các hệ số theo t = (x-{x0})/{h} (lũy thừa tăng dần) là:',
        file=output,
    )
    print(t_coeffs, file=output)

    print(
        f'Thay t = {t} vào đa thức theo t, ta được giá trị tại {x}:',
        file=output,
    )
    print(sum(i * t**m for m, i in enumerate(t_coeffs)), file=output)

    error = abs(ys[0] / fact * t * prod((t_sq - i**2 for i in range(1, n + 1))))
    print(file=output)
    print('Sai số ước lượng không vượt quá:', file=output)
    print(error, file=output)
    print('=' * 60, file=output)
def main():
    if not xs[0] <= start <= xs[-1]:
        exit(colored('Không hỗ trợ ngoại suy','red'))
    with open('output.txt', 'w', encoding='utf-8') as output:
        Stirling_overhaul(xs,ys,start,output)
    print('Đã ghi lời giải chi tiết vào file output.txt')
main()