from sympy import *
from math import *
import sys


class rungekutta_oop:
#{
    def __init__(self, expr, x_0, y_0, h, n, s):
    #{
        x          = symbols("x")
        y          = symbols("y")
        func       = sympify(expr)


        self.y_0   = y_0           # Giá trị ban đầu của x
        self.x_0   = x_0           # Giá trị ban đầu của y
        self.h     = h             # Giá trị h - khoảng cách giữa 2 điểm
        self.n     = n             # Số mốc cần đưa ra
        self.s     = s             # Số nấc s của phương pháp


        self.f     = lambdify((x,y), func, "math");
        self.sym_f = func;


        self.alpha = [0];
        self.beta  = [0];
        self.r     = [0];
    #}

    def __computeNext(self, prev):
    #{
        s     = 3 if(self.s == "heun" or self.s == "Heun") else self.s; 
        f     = self.f;
        h     = self.h;
        r     = self.r;
        alpha = self.alpha;
        beta  = self.beta;

        x   = prev[0]; 
        y   = prev[1]; 
        k   = [0];
        ret = y;

        for i in range(1, s+1):
        #{
            temp_x = x + alpha[i] * h;
            temp_y = y;

            for j in range(1, i):
            #{
                temp_y += beta[i-1][j] * k[j];
            #}
            
            k.append(h * f(temp_x, temp_y));
            ret  += r[i] * k[i];
        #}

        return ret;
    #}
    def __loadCoeff(self, s):
    #{
        if(s == 1):
        #{
            self.r     = [0, 1];
            self.alpha = [0, 0];
            self.beta  = [
                [0],
                [0, 0]
            ];
        #}
        elif(s == 2):
        #{
            self.r     = [0, 1/3, 2/3];
            self.alpha = [0, 0, 3/4];
            self.beta  = [
                [0],
                [0, 3/4]
            ];
        #}
        elif(s == 3):
        #{
            self.r     = [0, 1/6, 2/3, 1/6];
            self.alpha = [0, 0, 1/2, 1];
            self.beta  = [
                [0],
                [0, 1/2],
                [0, -1, 2]
            ];
        #}
        elif(self.s == "heun" or self.s == "Heun"):
        #{
            self.r     = [0, 1/4, 0, 3/4];
            self.alpha = [0, 0, 1/3, 2/3];
            self.beta  = [
                [0],
                [0, 1/3],
                [0, 0, 2/3]
            ];
        #}
        elif(s == 4):
        #{
            self.r     = [0, 1/6, 1/3, 1/3, 1/6];
            self.alpha = [0, 0, 1/2, 1/2, 1];
            self.beta  = [
                [0],
                [0, 1/2],
                [0, 0, 1/2],
                [0, 0, 0, 1]
            ];
        #}
    #}


    def Solve(self):
    #{
        self.__loadCoeff(self.s);
        grid = [[self.x_0, self.y_0]];

        for i in range(self.n):
        #{
            new_x = grid[-1][0] + self.h;
            new_y = self.__computeNext(grid[-1]);

            grid.append([new_x, new_y]);
        #}

        return grid;
    #}
#}







# Cấp chính xác
precision_x   = 3;
precision_y   = 7;
precision_eps = 12;


# Nhập dữ liệu
#y' = f(x,y)
expr      = "0.25*(x+1)*y*ln(50/y) "
x_0       = 0
y_0       = 2
h         = 0.1
n         = 4
# input("Nhập số nấc Runge-Kutta (1, 2, 3, 4 hoặc Heun):");
s         = 4

# Giải PTVP
s  = s if(s == "heun" or s == "Heun") else int(s); 
uu = rungekutta_oop(expr, x_0, y_0, h, n, s);
g  = (uu.Solve());



# In ra sai số và kết quả
if(s == "heun" or s == "Heun"):
#{
    print(f"Phương pháp Runge-Kutta 3 nấc với cấu hình hệ số kiểu Heun hoàn tất với lưới điểm sau, sai số toàn cục O(h^3) = {round(h**3, precision_eps)}:");
#}
else:
#{
    print(f"Phương pháp Runge-Kutta {s} nấc với hoàn tất với lưới điểm sau, sai số toàn cục O(h^{s}) = {round(h**s, precision_eps)}:");
#}
for x in g: print(round(x[0], precision_x), round(x[1], precision_y));



# Nghiệm kiểm chứng
print("====================================================================")
sol_expr  = input("Nhập phương trình nghiệm kiểm chứng y(x) hoặc NONE nếu không có: y(x) = ");
if(sol_expr != "NONE"):
#{
    sol = lambdify(symbols("x"), sympify(sol_expr), "math");
    print("So sánh với giá trị đúng của nghiệm:");
    for x in g: print(f"y*({round(x[0], precision_x)}) = {round(x[1], precision_y)}, so sánh với giá trị chuẩn y({round(x[0], precision_x)}) = {round(sol(x[0]), precision_y)}, sai số {round(abs(sol(x[0]) - x[1]), precision_eps)}");
#}



import matplotlib.pyplot as plt
def PlotPairs(pairList):
    t,x = zip(*pairList)
    plt.scatter(t,x)
    plt.show()

PlotPairs(g)