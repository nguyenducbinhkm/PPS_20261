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
def vietaFormula(coeff:list,x0:float) -> list: # nhận vào hệ số bậc tăng dần coeff và số thực x0, trả về list hệ số của tích 2 đa thức khi nhân vào
                                                #VD: vietaFormula([2,-3,1],3) -> [-6,11,-6,1]
    coeff = [0]+coeff
    n = len(coeff)
    for i in range(n-1):
        coeff[i]-=x0*coeff[i+1]
    return coeff # 3(n-1) mul
def plus_poly(high: list,low: list)->list: # cộng 2 đa thức, high là đa thức hệ số bậc cao hơn, low là đa thức hệ số bậc thấp hơn, trả về list hệ số
                                            #VD: plus_poly([1,2,3,4],[2,4,3]) -> [3,6,6,4]
    n = len(high)
    if n==1:
        return high
    for i in range(n-1):
        high[i]+=low[i]
    return high # len(high) cal
def mul_list(kirito: list,ratio: float) -> list: # VD: mul_list([1,2,3,4],1111) -> [1111,2222,3333,4444]
    return [i*ratio for i in kirito] # len(kirito) cal


def deri(ys:list): # trả về sai phân cấp tiếp theo của list ys
                    # VD: deri([1,5,2,3]) -> [4,-3,1]
    n = len(ys)
    return [ys[i+1]-ys[i] for i in range(n-1)]    # len(ys)*2 cal
def nodes(xs:list,c:float): #nhận vào list xs và số thực c, trả về
    p=len(xs)
    h=xs[1]-xs[0]
    j=0
    if abs(c-xs[0]) < abs(c-xs[-1]):
        i=0
        while abs(xs[i]-c)>=h:
            i+=1
        if xs[i]==c:
            exit(colored(f'VALUE KNOWN AT {c} IS: {ys[i]}','blue'))
        while xs[i+j+1]-xs[i-j]<=1 and 0<=i+j+1<=p-1 and 0<=i-j<=p-1:
            j+=1
        j-=1
        if j<=1:
            if (i-2)*(i+3)>=0:
                
                return range(i-2,i+4)
            else:
                exit(colored('Too few usable data points to perform Bessel Interpolation'.upper(),'red'))
        elif 2<=j<=3:
            
            return range(i-j,i+j+2)
        else:
            
            return range(i-3,i+5)
    else:
        i=-1
        if xs[i]==c:
            exit(colored(f'VALUE KNOWN AT {c} IS: {ys[i]}','blue'))
        while abs(xs[i]-c)>h:
            i-=1
        i-=1
        if xs[i]==c:
            exit(colored(f'VALUE KNOWN AT {c} IS: {ys[i]}','blue'))
        while xs[i+j+1]-xs[i-j]<=1 and -p<=i+j+1<=-1 and -p<=i-j<=-1:
            j+=1
        j-=1
        if j<=1:
            if (i-2)*(i+3)>0:
                
                return range(i-2,i+4)
            else:
                exit(colored('Too few usable data points to perform Bessel Interpolation'.upper(),'red'))
        elif 2<=j<=3:
            
            return range(i-j,i+j+2)
        else:
            
            return range(i-3,i+5)
def Bessel_overhaul(xs:list,ys:list,x:float,output):
    # NOTE: all the below listed calculations only consider the highest power, which is 2, and its coefficients

    ys = [ys[i] for i in nodes(xs,x)] ### phải cắt ys trước nếu không sẽ ảnh hưởng algorithm
    xs = [xs[i] for i in nodes(xs,x)]
    print(f'NỘI SUY BESSEL TẠI x = {x}', file=output)
    print('=' * 60, file=output)
    print('Các mốc được chọn và giá trị tương ứng:', file=output)
    for m,i in enumerate(xs):
        print(f'x = {xs[m]}, y = {ys[m]}', file=output)
    h = xs[1]-xs[0]
    n = int(len(xs)/2-1)
    x0=xs[n]
    y0 = ys[n]
    t = (x-xs[n])/h-1/2
    print(file=output)
    print(f'Bước h = {h}', file=output)
    print(f'Mốc trung tâm x0 = {x0}', file=output)
    print(f'Giá trị trung tâm y0 = {y0}', file=output)
    print(f't = (x - x0) / h - 1/2 = {t}', file=output)
    print(file=output)
    print('Các sai phân được sử dụng:', file=output)
    coeff=[-0.25,1]
    t_coeffs = []
    even_coeff = [(ys[n]+ys[n+1])/2]
    ys = deri(ys[:])
    odd_coeff = [ys[n]]
    print(f'Giá trị sai phân bậc 1 trung tâm = {odd_coeff[0]}', file=output)
    fact = 1
    for i in range(1,n+1):
        ys = deri(ys[:]) #n^2 cal total
        
        fact *= 2*i
        diff = (ys[n-i]+ys[n-i+1])/(2*fact)
        even_coeff = plus_poly(mul_list(coeff,diff),even_coeff) # n^2 cal
        print(f'Hệ số chẵn b{i} = {diff}', file=output)
        ys = deri(ys[:]) #n^2 cal
        
        fact *= 2*i+1
        diff = ys[n-i]/fact
        odd_coeff = plus_poly(mul_list(coeff,diff),odd_coeff) # n^2 cal
        print(f'Hệ số lẻ a{i} = {diff}', file=output)
        coeff = vietaFormula(coeff,(i+1/2)**2) #3n^2/2 cal
    
    odd_coeff = np.array(odd_coeff)
    even_coeff = np.array(even_coeff)
    t_sq = t**2
    q = 1
    t_vec = [[1]]
    for i in range(n):
        q *= t_sq
        t_vec.append([q])
    t_vec = np.array(t_vec)

    result = even_coeff@t_vec + t*(odd_coeff@ t_vec)
    print(file=output)
    print(f'Giá trị nội suy Bessel tại x = {x}:', file=output)
    print(result, file=output)
    for j in range(n+1):
            t_coeffs.append(even_coeff[j])
            t_coeffs.append(odd_coeff[j])
    print(file=output)
    print(f'Các hệ số theo t = (x-{x0})/{h} - 1/2 (lũy thừa tăng dần) là:', file=output)
    print(t_coeffs, file=output)
    
    print(f'Thay t = {t} vào đa thức theo t, ta được giá trị tại {x}:', file=output)
    print(sum(i*t**m for m,i in enumerate(t_coeffs)), file=output) #neglected (exists for the sake of substantiation)

    error = abs(ys[0]/fact*prod((t_sq-(i+0.5)**2 for i in range(n+1))))
    print(file=output)
    print('Sai số ước lượng không vượt quá:', file=output)
    print(error, file=output)
    print('=' * 60, file=output)
def main():
    if not xs[0] <= start <= xs[-1]:
        exit(colored('Không hỗ trợ ngoại suy','red'))
    with open('output.txt', 'w', encoding='utf-8') as output:
        Bessel_overhaul(xs,ys,start,output)
    print('Đã ghi lời giải chi tiết vào file output.txt')
main()