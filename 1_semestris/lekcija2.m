%% 1
syms x y(x)
atr = dsolve(diff(y,x) + y * tan(x) == 1 / cos(x))
%% 2
clearvars
syms x y(x)
vien = diff(y, x, 2)  + 2 * diff(y, x) - 3 * y == 6 * exp(x) + x;
atr = dsolve(vien)
atr_simplified = simplify(atr)
%% 3
clearvars
syms x y(x)
vien = diff(y, x) + y * cos(x) == sin(x) * cos(x)
atr = dsolve(vien, y(0) == 1)
%% 4 
clearvars
syms x y(x)
vien = diff(y, x, 2) - 2 * diff(y, x) == exp(x) * (x^2 + x - 3);
Dy = diff(y, x)
atr = dsolve(vien, y(0) == 2, Dy(0) == 2)
%% 5
clearvars
syms x y(x)
atr=dsolve((1+x)*diff(y, x))
%% patst 1. a
clearvars
syms x y(x)
vien = diff(y, x) + x * y ^ (1/3) == 3 * y;
atr = dsolve(vien)
%% patst 1. b
clearvars
syms x y(x)
vien = diff(y, x) - 2 * x * y == 3 * x ^ 3 * y ^ 2;
atr = dsolve(vien)
%% patst 1. c - kkas nesakr'it
clearvars
syms x y(x)
vien = sqrt(1 + y ^ 2) == -2 * (x^2 * y + y) * diff(y, x);
atr = dsolve(vien)
%% patst 2. a
clearvars
syms x y(x)
vien = diff(y, x, 3) == x * log(x);
atr = simplify(dsolve(vien))
%% patst 2. b
clearvars
syms x y(x)
vien = diff(y, x, 3) - diff(y, x, 2) - 2 * diff(y, x) == 4 * x + 3 * sin(x) + cos(x)
atr = simplify(dsolve(vien))
%% patst 2. c
clearvars
syms x y(x)
vien = diff(y, x, 2) + 2 * diff(y, x) + y == 3 * exp(-x)*sqrt(x + 1);
atr = simplify(dsolve(vien))
%% patst 3. a
clearvars
syms x y(x)
vien = diff(y, x)*sqrt(1 - x^2) + y == asin(x);
atr = dsolve(vien, y(0) == 0)
%% patst 3. b
clearvars
syms x y(x)
vien = diff(y, x) * cot(x) + y == 2 ;
atr = dsolve(vien, y(0) == -1)
%% patst 3. c
clearvars
syms x y(x)
vien = x * diff(y, x) + y == y ^ 2 * log(x);
atr = dsolve(vien, y(1) == 1)
%% patst 4. a
clearvars
syms x y(x)
vien = diff(y, x, 2) - 2 * diff(y,x) == exp(x) * (x ^ 2 + x - 3);
Dy = diff(y, x)
atr = simplify(dsolve(vien, y(1) == 2, Dy(0) == 2))
%% patst 4. b
clearvars
syms x y(x)
vien = diff(y, x, 2) + 6 * diff(y,x) + 9 *y == 10 * sin(x);
Dy = diff(y, x)
atr = dsolve(vien, y(0) == 0, Dy(0) == 0)
%% patst 4. c
clearvars
syms x y(x)
vien = diff(y, x, 2) + 6 * diff(y,x) + 9 *y == 10 * sin(x);
Dy = diff(y, x)
atr = dsolve(vien, y(0) == 0, Dy(0) == 0)